#include <iostream>
#include <vector>
#include <cmath>
#include <numbers>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <string>
using namespace std;
int main()
{
    int n = 1000000000;
    vector<int> lp(n + 1);
    vector<int> phi(n + 1);
    vector<int> pn;
    phi[1] = 1;
    for (int i = 2; i <= n; i++)
    {
        if (lp[i] == 0)
        {
            lp[i] = i;
            pn.push_back(i);
            phi[i] = i - 1;
        }

        for (int p : pn)
        {
            if (i * p > n)
                break;
            int newn = i * p;
            lp[newn] = p;
            if (i % p == 0)
            {
                phi[newn] = phi[i] * p;
                break;
            }
            else
            {
                phi[newn] = phi[i] * (p - 1);
            }
        }
    }
    vector<vector<long double>> allconst(1001, vector<long double>(1000001));

    int counter = 1;
    for (int i = 1; i <= n; i++)
    {
        counter = (i - 1) / 1000000 + 1;
        if (i < 3 && counter == 1)
        {
            allconst[counter][i % 1000000] = -1e100L;
            continue;
        }
        int locali = (i - 1) % 1000000;
        allconst[counter][locali] = ((long double)i / phi[i] - exp(numbers::egamma_v<long double>) * log(log(i))) * log(log(i));
    }
    vector<long double> maxconstonforeverydiapozons(1001, -1e100L);
    for (int i = 1; i <= 1000; i++)
    {
        maxconstonforeverydiapozons[i] = *max_element(allconst[i].begin(), allconst[i].end());
    }
    if (!filesystem::exists("graphs"))
    {
        filesystem::create_directory("graphs");
    }

    ofstream txt("graphs/maxima_report.txt");
    txt << "=== МАКСИМУМЫ allconst ПО МИЛЛИОННЫМ ИНТЕРВАЛАМ ===\n\n";
    txt << "Интервал\tМаксимум\tДостигнут на числе\n";

    ofstream csv("graphs/maxima_report.csv");
    csv << fixed << setprecision(13);
    csv << "Интервал;Максимум;Достигнут_на_числе\n";

    for (int block = 1; block <= 1000; block++)
    {
        int startNum = (block - 1) * 1000000 + 1;
        int endNum = block * 1000000;
        if (startNum < 3)
            startNum = 3;
        if (endNum > n)
            endNum = n;

        long double maxVal = -1e100L;
        int maxIndex = -1;

        for (int i = startNum; i <= endNum; i++)
        {
            int locali = (i - 1) % 1000000;
            if (allconst[block][locali] > maxVal)
            {
                maxVal = allconst[block][locali];
                maxIndex = i;
            }
        }

        string factorstr = "";
        if (maxIndex > 1)
        {
            int temp = maxIndex;
            while (temp > 1)
            {
                factorstr += to_string(lp[temp]) + " ";
                temp = temp / lp[temp];
            }
        }
        else
        {
            factorstr = "1";
        }

        txt << block << "-й миллион (" << startNum << "-" << endNum << ")\t"
            << maxVal << "\t"
            << maxIndex << "\t"
            << factorstr << "\n";

        csv << block << "-й миллион;"
            << maxVal << ";"
            << maxIndex << ";"
            << factorstr << "\n";
    }

    txt.close();
    csv.close();
    cout << "Created: graphs/maxima_report.txt and graphs/maxima_report.csv" << endl;

    ofstream html("graphs/maxima_per_million.html");
    html << "<!DOCTYPE html>\n<html>\n<head>\n"
         << "<meta charset='UTF-8'>\n"
         << "<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>\n"
         << "</head>\n<body style='background-color: #f5f5f5; font-family: sans-serif;'>\n";

    html << "<div id='chart' style='width:100%; max-width:1400px; height:800px; margin: 0 auto;'></div>\n<script>\n";

    html << "var x_vals = [";
    for (int block = 1; block <= 1000; block++)
    {
        html << block * 1000000 << (block == 1000 ? "" : ",");
    }
    html << "];\n";

    html << "var y_vals = [";
    for (int block = 1; block <= 1000; block++)
    {
        html << maxconstonforeverydiapozons[block] << (block == 1000 ? "" : ",");
    }
    html << "];\n";

    html << R"(
    var trace1 = {
        x: x_vals,
        y: y_vals,
        type: 'scattergl',
        mode: 'lines+markers',
        name: 'Максимум allconst',
        line: {color: 'blue', width: 2},
        marker: {color: 'red', size: 6}
    };

    var data = [trace1];

    var layout = {
        title: 'Максимальные значения allconst по миллионным интервалам (до 1 млрд)',
        autosize: true,
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#f5f5f5',
        showlegend: true,
        xaxis: { title: 'n (миллионы)', gridcolor: '#E1E1E1' },
        yaxis: { title: 'max(allconst)', gridcolor: '#E1E1E1' }
    };

    Plotly.newPlot('chart', data, layout);
    )" << "\n</script>\n</body>\n</html>";

    html.close();
    cout << "Created: graphs/maxima_per_million.html" << endl;
}