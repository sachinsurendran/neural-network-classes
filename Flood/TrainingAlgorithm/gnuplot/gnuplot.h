#ifndef Gnuplot_class
#define Gnuplot_class

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

class GnuPlot
{
    private:
        string gnuplot_data_file;
        fstream fp;
    public:
        int write(int x, double y);
        int write(int x, double y, double z);
	int write(int w, int x, double y, double z);
        GnuPlot();
        ~GnuPlot();

};

GnuPlot::GnuPlot()
{
    time_t rawtime;
    struct tm * timeinfo;
    char filename[100];

    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    strftime (filename, 100, "graph/%m_%d_%H_%M_plot.dat", timeinfo);
    gnuplot_data_file = filename;
    cout <<" Opening Gnuplot data file" << gnuplot_data_file << endl;
    fp.open(gnuplot_data_file.c_str(), fstream::in | fstream::out | fstream::app);
}

GnuPlot::~GnuPlot()
{
    fp.close();
}


int GnuPlot::write(int x, double y)
{
    cout << "Writing " << x <<":" << y  << endl;

    fp << x << " " << y << endl;
}

int GnuPlot::write(int x, double y, double z)
{
    cout << "Writing " << x << ":" << y << ":" << z << endl;
    fp <<x << " " << y << " " << z << endl;
}

int GnuPlot::write(int w, int x, double y, double z)
{
    cout << "Writing " << w << " : " << x << " : " << y << " : " << z << endl;
    fp << w << " " << x << " " << y << " " << z << endl;
}

#endif
