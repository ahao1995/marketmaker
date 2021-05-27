#include <deque>
#include <fstream>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include "IWCStrategy.h"
#include "Timer.h"
#include "as_utils.h"
#include "AverageVolatility.hpp"

USING_WC_NAMESPACE

struct order
{
    order() {}
    order(int id, double p, bool i, long time, LfOrderStatusType sts) : rid(id),
                                                                        price(p),
                                                                        is_buy(i),
                                                                        close_time(time),
                                                                        order_status(sts),
                                                                        cancel_cnt(0)
    {
    }
    int rid;
    long close_time;
    LfOrderStatusType order_status;
    int cancel_cnt; //撤销的次数
    double price;
    bool is_buy;
    bool operator<(const order &b) const
    {
        return close_time > b.close_time;
    }
};

class AsMakerStrategy : public IWCStrategy
{
public:
    virtual void init();
    virtual void on_market_tick(const LFTickMarketDataField *data, short source, long rcv_time);
    virtual void on_market_Dp5(const LFDepth5MarketDataField *data, short source, long rcv_time);
    virtual void on_market_trade(const LFTradeDataField *data, short source, long rcv_time);

    virtual void on_rsp_position(const PosHandlerPtr posMap, int request_id, short source, long rcv_time);
    virtual void on_rtn_trade(const LFCryptoRtnOrderField *data, int request_id, short source, long rcv_time);
    virtual void on_rtn_order(const LFCryptoRtnOrderField *data, int request_id, short source, long rcv_time);
    virtual void on_rsp_order(const LFCryptoInputOrderField *data, int request_id, short source, long rcv_time, int errorId = 0, const char *errorMsg = nullptr);

public:
    AsMakerStrategy(const string &name);
    ~AsMakerStrategy();
    //计算波动
    double get_volatility(double spread);
    //收集数据
    void collect_market_variables(double price, long rcv_time);
    //判断算法是否准备好
    bool is_algorithm_ready();
    //volatility 超过阈值 需要调整参数
    double volatility_diff_from_last_parameter_calculation(double cur_vol);
    //c_recalculate_parameters
    void recalculate_parameters(double price, double vol);
    // 计算报价
    void calculate_reserved_price_and_optimal_spread(double price, double vol);
    //处理过时的订单
    friend void cancel_aged_orders(AsMakerStrategy *);
    //计算剩余的钱
    double calculate_left_banlance();
    //当前价格是否有挂单
    bool has_qutoe(double ask_price, double bid_price, double cur_price);

    std::string instrument = "UNIUSDT"; //交易对
    double min_tick_price = 0.001;
    double trade_size = 1; //每次交易的数量
    double gamma = 0;
    // double sigma = 0.025;
    double kappa = 0;
    double m_eta;

    bool td_connected;   //通道是否就绪
    double init_balance; //剩余的钱

    int count = 0;

    std::unordered_map<int, order> m_orders;  //rid:对应的订单方便查看状态
    std::priority_queue<order> m_sort_orders; //根据时间进行排序的订单 先过期的订单在上面

private:
private:
    const double ONE_DAY_IN_MS = 86400000;
    const double ONE_HOUR_IN_MS = 3600000;
    double m_time_left;               // 时间窗口剩余
    double m_last_timestamp;          //最后一次更新时间
    double m_optimal_spread;          //计算得到的差
    double m_reserved_price;          //计算得到的中间价格r
    double m_optimal_ask;             //计算得到的ask报价
    double m_optimal_bid;             //计算得到的bid报价
    double m_froze_balance;           //被冻结的资金
    double m_inventory_risk_aversion; //

    AverageVolatility m_avg_vol;
    double m_latest_parameter_calculation_vol; //最后一次计算的vol
    bool should_qutoe;                         //应该发单

    /*需要传入的参数*/
    double m_closing_time;             //时间窗口大小
    double m_vol_to_spread_multiplier; //vol阈值
    bool m_parameters_based_on_spread; //价差变大要重新调节参数
    double m_min_spread;               //最小差百分比
    double m_max_spread;               //最大差百分比
    long m_aged_time;                  //订单过期时间
    double m_per_cnt;                  //合约每一张代表的面值 0.01uni一张
    int m_leverage;                    //杠杆倍数
    int m_delay_tick;                  //隔多少个tick发一次单
    int m_tick_cnt;                    //tick计数
    double m_order_refresh_tolerance_pct;
};

AsMakerStrategy::AsMakerStrategy(const string &name) : IWCStrategy(name),
                                                       m_avg_vol(180, 1),
                                                       m_closing_time(3600000),
                                                       m_min_spread(0.00005),
                                                       m_max_spread(0.01),
                                                       m_aged_time(10 * 1e9),
                                                       m_delay_tick(5),
                                                       m_order_refresh_tolerance_pct(0.05)
{
}

AsMakerStrategy::~AsMakerStrategy()
{
    //撤销全部订单
    util->cancel_all_order(SOURCE_BINANCE, "UNIUSDT");
}
void AsMakerStrategy::init()
{
    vector<string> ticker_binance;
    data->add_market_data(SOURCE_BINANCE);
    data->add_register_td(SOURCE_BINANCE);
    ticker_binance.push_back("uniusdt@depth5@500ms");
    ticker_binance.push_back("uniusdt@ticker");
    ticker_binance.push_back("uniusdt@aggTrade");
    util->subscribeMarketData(ticker_binance, SOURCE_BINANCE);

    //strategy related
    m_time_left = this->ONE_HOUR_IN_MS;
    init_balance = 0.1;
    m_per_cnt = 1;
    // m_per_cnt = 10;

    m_leverage = 5;
    m_froze_balance = 0;
    m_inventory_risk_aversion = 0.99;
    m_tick_cnt = 0;
    m_parameters_based_on_spread = true;
    m_vol_to_spread_multiplier = 1.3;
    m_latest_parameter_calculation_vol = 0;
}

void AsMakerStrategy::collect_market_variables(double price, long rcv_time)
{
    m_time_left = std::max(m_time_left - (rcv_time - m_last_timestamp) / 1000000, 0.0);
    m_avg_vol.add_sample(price);
    if (m_time_left == 0)
    {
        m_time_left = this->ONE_HOUR_IN_MS;
    }
}

double AsMakerStrategy::get_volatility(double spread)
{
    double vol = m_avg_vol.current_value();
    //TODO:细节补充
    if (eq(vol, 0))
    {
        if (!eq(m_latest_parameter_calculation_vol, 0))
        {
            vol = m_latest_parameter_calculation_vol;
        }
        else
        {
            vol = spread / 2;
        }
    }
    return vol;
}
bool AsMakerStrategy::is_algorithm_ready()
{
    return m_avg_vol.is_sampling_buffer_full();
}
void AsMakerStrategy::recalculate_parameters(double price, double vol)
{
    //计算gamma
    PosHandlerPtr pos = this->data->get_pos(SOURCE_BINANCE);
    pos->get_holding_net_pos("UNIUSDT");
    // position
    double q = pos->get_holding_net_pos("UNIUSDT");
    double min_spread = m_min_spread * price;
    double max_spread = m_max_spread * price;
    // double max_possible_gamma = std::min(
    //     (max_spread - min_spread) / (2 * std::abs(q) * (vol * vol)),
    //     (max_spread * (2 - m_inventory_risk_aversion) / m_inventory_risk_aversion + min_spread) / (vol * vol));
    // std::cout << "max_possible_gamma:" << max_possible_gamma << std::endl;
    // gamma = m_inventory_risk_aversion * max_possible_gamma;
    // if (gamma >= 0.001)
    //     gamma = 0.001;
    gamma = 1e-4;
    //计算k
    double max_spread_around_reserved_price = max_spread * (2 - m_inventory_risk_aversion) + min_spread * m_inventory_risk_aversion;
    std::cout << "max_spread_around_reserved_price" << max_spread_around_reserved_price << std::endl;
    if (max_spread_around_reserved_price <= gamma * (vol * vol))
    {
        kappa = std::exp(100);
    }
    else
    {
        kappa = gamma / (std::exp((max_spread_around_reserved_price * gamma - (vol * gamma) * (vol * gamma)) / 2) - 1);
    }
    printf("gamma:%.4lf,kappa:%4lf", gamma, kappa);
    //计算eta
    double q_where_to_decay_order_amount = q * (1 - m_inventory_risk_aversion);
    m_eta = 1;
    if (q_where_to_decay_order_amount != 0)
    {
        m_eta = m_eta / q_where_to_decay_order_amount;
    }
    m_latest_parameter_calculation_vol = vol;
}

double AsMakerStrategy::volatility_diff_from_last_parameter_calculation(double cur_vol)
{
    if (eq(m_latest_parameter_calculation_vol, 0))
    {
        return 0;
    }
    return std::abs(m_latest_parameter_calculation_vol - cur_vol) / m_latest_parameter_calculation_vol;
}

double AsMakerStrategy::calculate_left_banlance()
{
    PosHandlerPtr pos = this->data->get_pos(SOURCE_BINANCE);
    return init_balance * m_leverage * 0.95 - m_froze_balance;
}

void AsMakerStrategy::on_market_tick(const LFTickMarketDataField *md, short source, long rcv_time)
{
    // std::cout << "[tick] "
    //   << "InstrumentID:" << md->InstrumentID << " bid1:" << md->BidPrice1 << " ask1:" << md->OfferPrice1 << std::endl;
}

void AsMakerStrategy::calculate_reserved_price_and_optimal_spread(double price, double vol)
{
    double time_left_fraction = m_time_left / m_closing_time;
    PosHandlerPtr pos = this->data->get_pos(SOURCE_BINANCE);
    pos->get_holding_net_pos("UNIUSDT");
    // position
    double q = pos->get_holding_net_pos("UNIUSDT");
    double mid_price_variance = std::pow(vol, 2);

    m_reserved_price = price - q * this->gamma * mid_price_variance * time_left_fraction;
    m_optimal_spread = this->gamma * mid_price_variance * time_left_fraction + 2 * std::log(1 + this->gamma / this->kappa) / this->gamma;

    //限制价差在合理范围内，价差要大于最小价差 小于最大价差
    double min_limit_bid = std::min(price * (1 - m_max_spread), price - m_vol_to_spread_multiplier * vol);
    double max_limit_bid = price * (1 - m_min_spread);
    double min_limit_ask = price * (1 + m_min_spread);
    double max_limit_ask = std::max(price * (1 + m_max_spread), price + m_vol_to_spread_multiplier * vol);

    m_optimal_ask = std::min(std::max(m_reserved_price + m_optimal_spread / 2, min_limit_ask),
                             max_limit_ask);
    m_optimal_bid = std::min(std::max(m_reserved_price - m_optimal_spread / 2, min_limit_bid),
                             max_limit_bid);
}

void cancel_aged_orders(AsMakerStrategy *str)
{
    long cur_time = kungfu::yijinjing::getNanoTime();
    while (!str->m_sort_orders.empty())
    {
        auto cur_order = str->m_sort_orders.top();
        auto &cur_sts = str->m_orders[cur_order.rid].order_status;
        if (cur_order.close_time <= cur_time)
        {
            //已经取消了的订单、已经过期、已经成交的订单直接erase
            if (cur_sts == LF_CHAR_Canceled || cur_sts == LF_CHAR_OrderExpired || cur_sts == LF_CHAR_AllTraded)
            {
                // std::cout << "pop:" << cur_order.rid << std::endl;
                str->m_orders.erase(cur_order.rid);
                str->m_sort_orders.pop();
            }
            //未成交的撤单
            else if (cur_sts == LF_CHAR_OrderInserted)
            {
                str->util->cancel_order(SOURCE_BINANCE, cur_order.rid);
                cur_sts = LF_CHAR_OrderCancelSent;
                str->m_orders[cur_order.rid].cancel_cnt += 1;
                str->m_sort_orders.pop();
                //插入一个新订单
                str->m_sort_orders.push(order(cur_order.rid, cur_order.price, cur_order.is_buy, cur_time + 15 * 1e9, LF_CHAR_OrderCancelSent));
            }
            //对于第二次撤单的来说
            else if (cur_sts == LF_CHAR_OrderCancelSent)
            {
                //撤单次数小于3次就继续撤单，否则认为该订单已经被撤掉了
                if (str->m_orders[cur_order.rid].cancel_cnt <= 5)
                {
                    str->util->cancel_order(SOURCE_BINANCE, cur_order.rid);
                    str->m_orders[cur_order.rid].cancel_cnt += 1;
                    str->m_sort_orders.pop();
                    //插入一个新订单
                    str->m_sort_orders.push(order(cur_order.rid, cur_order.price, cur_order.is_buy, cur_time + 15 * 1e9, LF_CHAR_OrderCancelSent));
                }
                else
                {
                    str->m_orders.erase(cur_order.rid);
                    str->m_sort_orders.pop();
                }
            }
            else if (cur_sts == LF_CHAR_OrderSent)
            {
                str->m_sort_orders.pop();
                //插入一个新订单
                str->m_sort_orders.push(order(cur_order.rid, cur_order.price, cur_order.is_buy, cur_time + 3 * 1e9, LF_CHAR_OrderSent));
            }
        }
        else
        {
            break;
        }
    }
}

bool AsMakerStrategy::has_qutoe(double ask_price, double bid_price, double cur_price)
{

    bool flag = false;
    for (auto &cur_order : m_orders)
    {
        if (cur_order.second.is_buy == true &&
            eq(cur_order.second.price, bid_price) &&
            cur_order.second.order_status != LF_CHAR_OrderCancelSent)
            flag = true;
        else if (cur_order.second.is_buy == false &&
                 eq(cur_order.second.price, ask_price) &&
                 cur_order.second.order_status != LF_CHAR_OrderCancelSent)
            flag = true;
        if (cur_order.second.order_status == LF_CHAR_OrderInserted &&
            (cur_order.second.price - cur_price) / cur_order.second.price > m_order_refresh_tolerance_pct)
        {
            util->cancel_order(SOURCE_BINANCE, cur_order.second.rid);
            cur_order.second.order_status = LF_CHAR_OrderCancelSent;
            cur_order.second.cancel_cnt += 1;
        }
    }
    return flag;
}
//TODO:
//挂单中含有了同样价格的那么就不挂单
//如果价格在撤单之前波动过大，那么直接撤掉不等待时间到达 每个tick检测一次（优化项）
//正在活动的订单的价格有个过期时间，30s以内不会再次以这个价格挂单
//用python 对日志文件进行分析 提取需要的数据-->csv

void AsMakerStrategy::on_market_Dp5(const LFDepth5MarketDataField *md, short source, long rcv_time)
{
    double mid_price = (md->BidPrice1 + md->OfferPrice1) / 2;
    collect_market_variables(mid_price, rcv_time);
    if (m_tick_cnt == m_delay_tick)
    {
        should_qutoe = true;
        m_tick_cnt = 0;
    }
    else
    {
        m_tick_cnt++;
    }
    if (is_algorithm_ready() && should_qutoe)
    {
        should_qutoe = false;
        double vol = get_volatility(std::abs(md->BidPrice1 - md->OfferPrice1));
        //是否需要重新计算参数
        if ((gamma - 0 < 1e-7 || kappa - 0 < 1e-7) ||
            (m_parameters_based_on_spread &&
             volatility_diff_from_last_parameter_calculation(vol) > (m_vol_to_spread_multiplier - 1)))
        {
            recalculate_parameters(mid_price, vol);
            KF_LOG_DEBUG_FMT(logger, "[recalculate_parameters] (gamma) %lf (kappa) %lf", gamma, kappa);
        }

        calculate_reserved_price_and_optimal_spread(mid_price, vol);
        //适应精度
        m_optimal_ask = round_nearest(m_optimal_ask, 0.0001);
        m_optimal_bid = round_nearest(m_optimal_bid, 0.0001);

        KF_LOG_DEBUG_FMT(logger, "[qutoe] (mid_price) %lf (ask1) %lf (bid1) %lf (reserved_price) %lf (optimal_spread) %lf (optimal_ask) %lf (optimal_bid) %lf (volatility) %lf (gamma) %lf (kappa) %lf",
                         mid_price,
                         md->OfferPrice1,
                         md->BidPrice1,
                         m_reserved_price,
                         m_optimal_spread,
                         m_optimal_ask,
                         m_optimal_bid,
                         vol,
                         gamma,
                         kappa);

        // printf("can_buy_uni: %.4lf\n", can_buy_uni);
        //判断价格是否合理 并且 开单有没有剩余 该价格有没有挂单
        double can_buy_uni = calculate_left_banlance() / m_optimal_ask;
        // printf("can_buy_uni: %.4lf\n", can_buy_uni);
        bool f = has_qutoe(m_optimal_ask, m_optimal_bid, mid_price);
        if (can_buy_uni > 2 * m_per_cnt && !f)
        {
            //下单 给订单设置个过期时间
            int rid_bid = util->insert_maker_order(SOURCE_BINANCE, "UNIUSDT", "binance",
                                                   m_optimal_bid, m_per_cnt,
                                                   LF_CHAR_Buy, LF_CHAR_Non);
            int rid_ask = util->insert_maker_order(SOURCE_BINANCE, "UNIUSDT", "binance",
                                                   m_optimal_ask, m_per_cnt,
                                                   LF_CHAR_Sell, LF_CHAR_Non);

            order order_bid(rid_bid, m_optimal_bid, true, rcv_time + m_aged_time, LF_CHAR_OrderSent);
            order order_ask(rid_ask, m_optimal_ask, false, rcv_time + m_aged_time, LF_CHAR_OrderSent);
            m_sort_orders.push(order_bid);
            m_sort_orders.push(order_ask);
            m_orders.insert(std::pair<int, order>(rid_bid, order_bid));
            m_orders.insert(std::pair<int, order>(rid_ask, order_ask));
            BLCallback bl = std::bind(cancel_aged_orders, this);
            util->insert_callback(kungfu::yijinjing::getNanoTime() + m_aged_time, bl);
        }
        else //持仓市价全平
        {
            //判断下所持的仓位是 多 还是 空

            // int rid_bid = util->insert_maker_order(SOURCE_BINANCE, "UNIUSDT", "binance",
            //                                        m_optimal_bid, m_per_cnt,
            //                                        LF_CHAR_Buy, LF_CHAR_Non);
            // int rid_ask = util->insert_maker_order(SOURCE_BINANCE, "UNIUSDT", "binance",
            //                                        m_optimal_ask, m_per_cnt,
            //                                        LF_CHAR_Sell, LF_CHAR_Non);
        }
    }
    m_last_timestamp = rcv_time;
}

void AsMakerStrategy::on_market_trade(const LFTradeDataField *md, short source, long rcv_time)
{
    // std::cout << "[trade] "
    //           << "InstrumentID:" << md->InstrumentID << " price:" << md->LastPrice << std::endl;
}

void AsMakerStrategy::on_rsp_position(const PosHandlerPtr posMap, int request_id, short source, long rcv_time)
{
    if (request_id == -1 && (source == SOURCE_BINANCE))
    {
        td_connected = true;
        KF_LOG_INFO(logger, "td connected");
        if (posMap.get() == nullptr)
        {
            data->set_pos(PosHandlerPtr(new PosHandler(source)), source);
        }
    }
    else
    {
        KF_LOG_DEBUG(logger, "[RSP_POS] " << posMap->to_string());
    }
}

void AsMakerStrategy::on_rsp_order(const LFCryptoInputOrderField *order, int request_id, short source, long rcv_time, int errorId, const char *errorMsg)
{
    if (errorId != 0)
        KF_LOG_ERROR(logger, " (err_id)" << errorId << " (err_msg)" << errorMsg << "(order_id)" << request_id << " (source)" << source);
}

//处理部分成交
void AsMakerStrategy::on_rtn_trade(const LFCryptoRtnOrderField *data, int request_id, short source, long rcv_time)
{
    m_orders[request_id].order_status = data->OrderStatus;

    KF_LOG_DEBUG_FMT(logger, "[on_rtn_trade] (id)%d (ref)%s (ticker)%s (Vsum)%lf (Vtrd)%lf (Status)%c,(avgPrice)%lf",
                     request_id,
                     data->OrderRef,
                     data->InstrumentID,
                     data->VolumeTotalOriginal,
                     data->VolumeTraded,
                     data->OrderStatus,
                     data->AvgPrice);
}

void AsMakerStrategy::on_rtn_order(const LFCryptoRtnOrderField *data, int request_id, short source, long rcv_time)
{
    m_orders[request_id].order_status = data->OrderStatus;
    if (data->OrderStatus == LF_CHAR_OrderInserted)
        m_froze_balance += data->LimitPrice * data->VolumeTotalOriginal;
    else if (data->OrderStatus == LF_CHAR_Canceled || data->OrderStatus == LF_CHAR_OrderExpired)
        m_froze_balance -= data->LimitPrice * data->VolumeTotalOriginal;
    KF_LOG_DEBUG_FMT(logger, "[on_rtn_order] (id)%d (ref)%s (ticker)%s (Vsum)%lf (Vtrd)%lf (Status)%c (limitPrice)%lf",
                     request_id,
                     data->OrderRef,
                     data->InstrumentID,
                     data->VolumeTotalOriginal,
                     data->VolumeTraded,
                     data->OrderStatus,
                     data->LimitPrice);
}

int main(int argc, const char *argv[])
{
    AsMakerStrategy str(string("as_maker"));
    str.init();
    str.start();
    str.block();
    return 0;
}
