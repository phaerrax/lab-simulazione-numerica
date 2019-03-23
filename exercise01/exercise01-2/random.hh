#ifndef __Random__
#define __Random__

class Random
{
    private:
        int m1, m2, m3, m4, l1, l2, l3, l4, n1, n2, n3, n4;

    public:
        Random();
        ~Random();
        //
        void SetRandom(int *, int, int);
        void SaveSeed();
        double Rannyu(void);
        double Rannyu(double min, double max);
        double Gauss(double mean, double sigma);
        double exponential(double rate);
        double cauchylorentz(double median, double scale);
};

#endif // __Random__
