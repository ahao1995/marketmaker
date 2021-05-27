#include <queue>
#include <cmath>
#include <iostream>

class AverageVolatility
{
public:
    AverageVolatility(int sampling_length, int processing_length)
    {
        m_sampling_length = sampling_length;
        m_processing_length = processing_length;
        m_sampling_sum = 0;
        m_processing_sum = 0;
    }
    void add_sample(double value)
    {
        if (!is_sampling_buffer_full())
        {
            m_sampling_sum += value;
            m_sampling_buffer.push_back(value);
        }
        else
        {
            m_sampling_sum -= m_sampling_buffer.front();
            m_sampling_buffer.pop_front();
            m_sampling_buffer.push_back(value);
            m_sampling_sum += value;
        }

        double indicator = indicator_calculation();
        if (!is_processing_buffer_full())
        {
            m_processing_sum += indicator;
            m_processing_buffer.push_back(indicator);
        }
        else
        {
            m_processing_sum -= m_processing_buffer.front();
            m_processing_buffer.pop_front();

            m_processing_buffer.push_back(indicator);
            m_processing_sum += indicator;
        }
    }
    double current_value()
    {
        return processing_calculation();
    }

    bool is_sampling_buffer_full()
    {
        return m_sampling_buffer.size() >= m_sampling_length;
    }

    bool is_processing_buffer_full()
    {
        return m_processing_buffer.size() >= m_processing_length;
    }

private:
    int m_sampling_length;
    std::deque<double> m_sampling_buffer;
    int m_processing_length;
    std::deque<double> m_processing_buffer;
    double m_sampling_sum;
    double m_processing_sum;

private:
    double indicator_calculation()
    {
        double mean = m_sampling_sum / m_sampling_buffer.size();
        double sum_temp = 0;
        if (m_sampling_buffer.size() == 1)
            return 0;
        for (int i = 0; i < m_sampling_buffer.size(); i++)
        {
            double last = m_sampling_buffer.back();
            sum_temp += std::pow(last - mean, 2);
            m_sampling_buffer.pop_back();
            m_sampling_buffer.push_front(last);
        }

        return (sum_temp / (m_sampling_buffer.size()));
    }

    double processing_calculation()
    {
        return std::sqrt(m_processing_sum / m_processing_buffer.size());
    }
};
