#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <iomanip>
#include <string>
#include <algorithm>   // для std::max, std::min
#include <limits>      // для std::numeric_limits

using namespace std;

const long long N = 4294967296LL;   // 2^32
const long long SEGMENT_SIZE = 1000000;
const long long L_START = 262;      // n0+1

const long double GAMMA = 0.5772156649015328606L;
const long double E_GAMMA = exp(GAMMA);

vector<int> sieve_primes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i * i <= limit; ++i) {
        if (is_prime[i]) {
            for (int j = i * i; j <= limit; j += i)
                is_prime[j] = false;
        }
    }
    vector<int> primes;
    for (int i = 2; i <= limit; ++i)
        if (is_prime[i]) primes.push_back(i);
    return primes;
}

// Функция для сохранения интерактивного HTML-графика
void savePlotHTML(const string& filename,
                  const vector<pair<double,double>>& points) {
    ofstream html(filename);
    if (!html) {
        cerr << "Cannot create " << filename << endl;
        return;
    }

    html << "<!DOCTYPE html>\n<html>\n<head>\n";
    html << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n";
    html << "</head>\n<body>\n";
    html << "<div id='chart' style='width:100%; max-width:1400px; height:900px; margin:0 auto;'></div>\n";
    html << "<script>\n";

    // x_vals
    html << "var x_vals = [";
    for (size_t i = 0; i < points.size(); ++i) {
        html << points[i].first;
        if (i != points.size() - 1) html << ",";
    }
    html << "];\n";

    // y_vals
    html << "var y_vals = [";
    for (size_t i = 0; i < points.size(); ++i) {
        html << setprecision(15) << points[i].second;
        if (i != points.size() - 1) html << ",";
    }
    html << "];\n";

    html << R"(
    var trace = {
        x: x_vals,
        y: y_vals,
        type: 'scattergl',
        mode: 'lines+markers',
        name: 'C_k (максимум f(n))',
        line: {color: 'blue', width: 1},
        marker: {color: 'red', size: 3}
    };

    var layout = {
        title: 'Кусочно-оптимальные константы для 2^32',
        xaxis: {title: 'Номер миллионного интервала'},
        yaxis: {title: 'Максимальное значение f(n)'},
        hovermode: 'closest'
    };

    Plotly.newPlot('chart', [trace], layout);
    )";

    html << "</script>\n</body>\n</html>\n";
}

int main() {
  
    ofstream table_file("table_c_k_full.csv");
    ofstream log_file("progress.log");
    if (!table_file.is_open() || !log_file.is_open()) {
        cerr << "Cannot open output files" << endl;
        return 1;
    }


    table_file << "interval_id;L;R;n;phi;f\n";
    table_file << setprecision(20);


    auto total_start = chrono::steady_clock::now();
    log_file << "Started full enumeration for N = " << N << endl;


    int sqrtN = static_cast<int>(sqrt(N)) + 1;
    vector<int> primes = sieve_primes(sqrtN);
    log_file << "Generated " << primes.size() << " primes up to " << sqrtN << endl;

    vector<pair<double,double>> plot_points;  // для графика

    long long L = L_START;
    int interval_id = 1;
    int last_progress_percent = -1;

    while (L <= N) {
        long long R = min(L + SEGMENT_SIZE - 1, N);
        size_t len = R - L + 1;

        // Массивы для сегмента
        vector<long long> phi(len), rem(len);
        for (size_t i = 0; i < len; ++i) {
            long long x = L + i;
            phi[i] = x;
            rem[i] = x;
        }


        for (int p : primes) {
            long long pp = p;
            if (pp * pp > R) break;
            long long first = ((L + pp - 1) / pp) * pp;
            if (first < pp * pp) first = pp * pp;
            for (long long x = first; x <= R; x += pp) {
                size_t idx = x - L;
                phi[idx] -= phi[idx] / p;
                while (rem[idx] % p == 0) rem[idx] /= p;
            }
        }
        for (size_t i = 0; i < len; ++i) {
            if (rem[i] > 1) {
                phi[i] -= phi[i] / rem[i];
            }
        }

        // Поиск максимума f(n) в сегменте
        long long local_best_n = L;
        long double local_best_f = -1e18L;
        for (size_t i = 0; i < len; ++i) {
            long long n = L + i;
            long double t = log(log((long double)n));
            long double r = (long double)n / (long double)phi[i];
            long double f = (r - E_GAMMA * t) * t;
            if (f > local_best_f) {
                local_best_f = f;
                local_best_n = n;
            }
        }


        table_file << interval_id << ";" << L << ";" << R << ";"
                   << local_best_n << ";" << phi[local_best_n - L] << ";"
                   << local_best_f << "\n";


        plot_points.emplace_back(interval_id, (long double)local_best_f);


        int progress = (int)(((long double)(R - L_START + 1) / (long double)(N - L_START + 1)) * 100);
        if (progress / 5 > last_progress_percent / 5) {
            last_progress_percent = progress;
            auto now = chrono::steady_clock::now();
            chrono::duration<double> elapsed = now - total_start;
            cout << setw(3) << progress << "% done; interval = " << interval_id
                 << "; time = " << elapsed.count() << " s" << endl;
            log_file << setw(3) << progress << "% done; interval = " << interval_id
                     << "; time = " << elapsed.count() << " s" << endl;
        }

        L += SEGMENT_SIZE;
        ++interval_id;
    }

    auto total_end = chrono::steady_clock::now();
    chrono::duration<double> total_elapsed = total_end - total_start;

    cout << "Finished in " << total_elapsed.count() << " seconds" << endl;
    log_file << "Finished in " << total_elapsed.count() << " seconds" << endl;

    table_file.close();
    log_file.close();


    savePlotHTML("plot.html", plot_points);
    cout << "Interactive plot saved to plot.html" << endl;

    return 0;
}
