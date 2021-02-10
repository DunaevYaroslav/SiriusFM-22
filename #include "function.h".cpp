#include "function.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <time.h>

#define AS(i, j) a[index[i] * n + (j)]
#define XS(i) x[index[i]]
#define BS(i) b[index[i]]
const double const1000 = 1000.0;
const double eps = 1e-50;
#define GIGA_MODIFIER 1e9
#define KILO_MODIFIER 1e3
#define MICRO_MODIFIER 1e-6
using namespace std;

double currentTimeNano()
{
    struct timespec t;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
    return (double)(t.tv_sec * GIGA_MODIFIER + t.tv_nsec);
}

class CCommonData
{
  public:
    double *a; // matrix pointer
    double *b; //
    int *index;
    int n; // total size of matrix
    int i; // current step
};
class CThreadData
{
  public:
    CCommonData *CD;
    pthread_t id; // Thread identifier
    int number;
    int rs;
    int j0;
    int ThreadCount;
    int re; // end row
    int n;
    double *x;
    double *a;
    double *b;
    double nevyazka;
    double sum;
    int count1;
};

void *Nev(void *Ptr)
{
    CThreadData *TD = (CThreadData *)Ptr;
    // int number = TD->number;
    int n = TD->n;
    int j0 = TD->j0;
    double *a = TD->CD->a; // matrix
    double *b = TD->CD->b;
    double *x = TD->x;
    int kk = TD->ThreadCount;
    int i;

    for (i = 0; i < n; i++) {
        if (i % kk == TD->number) {
            TD[i % kk].nevyazka += (a[j0 * n + i] * (x[i]) - b[i]) *
                                   (a[j0 * n + i] * (x[i]) - b[i]);
        }
    }

    return NULL;
}

void *Subtraction(void *Ptr)
{
    CThreadData *TD = (CThreadData *)Ptr;
    int number = TD->number;
    int n = TD->n;
    int j0 = TD->j0;
    double *a = TD->CD->a; // matrix
    double *b = TD->CD->b;
    int ThreadCount = TD->ThreadCount;
    int i;
    int j;
    for (j = j0 + 1; j < n; j++) {
        if (j % ThreadCount == number) {
            double norm = a[j * n + j0];
            for (i = j0; i < n; i++) {
                a[j * n + i] -= norm * a[j0 * n + i];
            }
            b[j] -= norm * b[j0];
        }
    }
    return NULL;
}

bool SolveSystem(int n, double *a, double *b, double *x, int kk, double &nevv)
{
    int i;
    int j;
    int k;
    int j0 = 0;
    double tmp;
    double *cop = new double[n * n];
    double norm;
    CCommonData CD;
    CD.a = a;
    CD.b = b;
    // CD.index = index;
    CD.n = n;
    // int i;
    int t;

    CThreadData *T = new CThreadData[kk];
    for (t = 0; t < kk; t++) {
        T[t].ThreadCount = kk;
        T[t].CD = &CD;
        T[t].number = t;
        T[t].n = n;
    }

    for (i = 0; i < n; i++) {
        x[i] = b[i];
        for (j = 0; j < n; j++) {
            cop[i * n + j] = a[i * n + j];
        }
    }
    double *time = new double[kk];
    for (i = 0; i < kk; i++) {
        time[i] = 0;
        T[i].count1 = 0;
    }
    double start1 = clock();
    for (j0 = 0; j0 < n; j0++) { // или же n-1
        for (i = j0; i < n; i++) {
            if (abs(a[i * n + j0]) > abs(a[j0 * n + j0])) { //замена
                for (k = j0; k < n; k++) {

                    tmp = a[i * n + k];
                    a[i * n + k] = a[j0 * n + k];
                    a[j0 * n + k] = tmp;
                }

                tmp = b[i];
                b[i] = b[j0];
                b[j0] = tmp;
            }
        }
        //нормализация
        norm = a[j0 * n + j0];

        if (abs(norm) < eps) {

            delete[] time;
            delete[] cop;
            delete[] T;
            return false;
        }

        for (k = j0; k < n; k++) {
            a[j0 * n + k] /= norm;
        }

        b[j0] /= norm;
        //вычитание параллельное
        int pp;
        int count = 0;
        for (pp = 0; pp < kk; pp++) {
            T[pp].j0 = j0;
        }

        T[pp].ThreadCount = kk;
        for (k = min(j0 + 1, n - 1); k < min(n, j0 + 1 + kk); k++) {
            count++;

            int f = k % kk;
            T[f].re = k;

            double all_time = currentTimeNano();
            pthread_create(&(T[f].id), 0, Subtraction, &T[f]);

            all_time = currentTimeNano() - all_time;
            time[f] += all_time;
        }
        for (int t = 0; t < kk; t++) { // t<count
            pthread_join(T[t].id, 0);
        }
    }

    for (i = 0; i < kk; i++) {
        cout << "Time of "
             << "thread " << i + 1 << ": "
             << time[i] / (const1000 * const1000 * const1000) << '\n';
    }
    if (n > 2) {
        for (int rr = n - 2; rr >= 0; rr--) {
            for (j = rr + 1; j < n; j++) {
                b[rr] = b[rr] - a[rr * n + j] * b[j];
            }
        }
    } else {
        b[1] = a[3];
        b[0] = a[0] - a[1] * a[3];
    }

    //вычислили решение
    //вычислим невязку
    double nevyazka = 0;
    for (int p = 0; p < kk; p++) {
        T[p].sum = 0;
        T[p].nevyazka = nevyazka;
    }
    for (int p = 0; p < kk; p++) {
        T[p].a = cop;
        T[p].b = x;
        T[p].x = b;
        for (i = 0; i < n; i++) {
            if (i % kk == p) {
                T[p].j0 = i;
                pthread_create(&(T[p].id), 0, Nev, &T[p]);
            }
        }
    }
    double finish1 = clock();
    for (int t = 0; t < kk; t++) {
        pthread_join(T[t].id, 0);
    }

    double nev = 0;
    for (i = 0; i < kk; i++) {
        nev += (T[i].nevyazka);
    }
    double nn = 0;
    for (i = 0; i < n; i++) {
        nn += (x[i]) * (x[i]);
    }
    nn = sqrt(nn);
    nev = sqrt(nev);

    nevv = nev / nn;
    delete[] time;
    delete[] T;
    delete[] cop;
    cout << "Time: " << (finish1 - start1) / (const1000 * const1000)
         << " sec\n";
    return true;
}


///////////////////////////////////////////////////////////////////////////

#include "function.h"
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <thread>
using namespace std;
int const NN = 1000;
int const const6 = 6;
int const const5 = 5;
int const const4 = 4;

double f(int k, int n, int i, int j)
{
    if (k == 1) {
        return double(n - max(i, j) + 1.0);
    }
    if (k == 2) {
        return double(max(i, j));
    }
    if (k == 3) {
        return double(abs(i - j));
    }
    if (k == 4) {
        return double(1.0 / (i + j - 1));
    }
    return 0;
}

double err(const double *b, int n)
{
    int i;
    double norm[n + 1];
    double s = 0;
    for (i = 0; i < n; i++) {
        norm[i] = (i + 1) % 2;
    }
    for (i = 0; i < n; i++) {
        s += (b[i] - norm[i]) * (b[i] - norm[i]);
    }
    return sqrt(s);
}

double nev(const double *x, const double *b, const double *a, int n)
{
    int i;
    int j;
    double Ax[n + 1];
    double sq1 = 0;
    double sq2 = 0;
    for (i = 0; i < n; i++) {
        Ax[i] = 0;
        for (j = 0; j < n; j++) {
            Ax[i] += a[i * n + j] * x[j];
        }
    }
    for (i = 0; i < n; i++) {
        sq1 += (Ax[i] - b[i]) * (Ax[i] - b[i]);
        sq2 += (b[i]) * (b[i]);
    }
    return (sqrt(sq1)) / (sqrt(sq2));
    //||Ax-b||/||b||
}

void create_b(int n, double *b, double **a)
{
    int i;
    int j;
    if (n % 2 == 0) {
        for (i = 1; i <= n; i++) {
            for (j = 1; j < n; j++) {
                if ((j % 2) == 1) {
                    b[i] += a[i][j];
                }
            }
        }
    } else {
        for (i = 1; i <= n; i++) {
            b[i] = 0;
            for (j = 1; j <= n; j++) {
                if (j % 2 == 1) {
                    b[i] += a[i][j];
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{

    int n = 0;
    int m = 0;
    int k = 0;
    int i;
    int j;
    int ss[const5];
    string filename;
    string s;

    if (argc < const4) {
        cout << "not enough numbers";
        return -1;
    }

    for (i = 1; i < const5; i++) {
        s = argv[i];
        ss[i] = s.size();
    }
    for (i = 1; i < const5; i++) {
        for (j = 0; j < ss[i]; j++) {
            if (isdigit(argv[i][j]) == 0) {
                cout << "not a number\n";
                return -1;
            }
        }
    }
    int kk = (atoi(argv[1]));
    n = (atoi(argv[2]));
    m = (atoi(argv[3]));
    k = (atoi(argv[4]));
    if (kk <= 0) {
        cout << "troubles with number of threads";
        return -1;
    }
    if (k > const4) {
        cout << "k very large\n";
        return -1;
    }
    if (n == 0) {
        cout << "matrix size = 0\n";
        return -1;
    }
    if (n < 0) {
        cout << "not a natural number\n";
        return -1;
    }
    if (m > n) {
        cout << "m > n\n";
        return -1;
    }
    if (n > NN) {
        cout << "very big n";
        return -2;
    }
    if (k == 0) {
        if (argc < const6) {
            cout << "file not defined";
            return -1;
        }
        filename = argv[const5];
    }
    double **a = new double *[n + 1];
    for (i = 0; i < n + 1; i++) {
        a[i] = new double[n + 1];
    }

    ifstream fin(filename);
    if (!fin.is_open() && (!(k != 0))) {
        cout << "Filename is bad\n";
        for (i = 0; i < n + 1; i++) {
            delete[] a[i];
        }
        delete[] a;
        return -3;
    }

    if (k == 0) {
        if (fin.is_open()) {
            int count = 0;
            int y = 0;
            for (int i = 1; i <= n; i++) {
                for (int j = 1; j <= n; j++) {
                    if (fin >> a[i][j]) {
                        if (!(fin.fail())) {
                            count++;
                        }
                    }
                    if (fin.fail() && count != 0) {
                        cout << "char in file";
                        for (y = 0; y < n + 1; y++) {
                            delete[] a[y];
                        }
                        delete[] a;
                        return -3;
                    }
                }
                if (count == 0) {
                    for (i = 0; i < n + 1; i++) {
                        delete[] a[i];
                    }
                    delete[] a;
                    cout << "empty file or only chars in file\n";
                    return -3;
                }
            }
        }
    }

    if (k != 0) {
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                a[i][j] = f(k, n, i, j);
            }
        }
    }

    double *b = new double[n + 1];
    double *copb = new double[n + 1];
    double *x = new double[n + 1];
    create_b(n, b, a);
    for (i = 1; i <= n; i++) {
        x[i] = b[i];
        copb[i] = b[i];
    }
    cout << 'A' << '\n';
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= m; j++) {
            printf("%.3e", a[i][j]);
            if (j != m) {
                cout << " ";
            }
        }
        printf(" %.3e", b[i]);
        cout << '\n';
        cout << " ";
    }

    srand(time(0));

    double *aa = new double[(n) * (n)];
    for (i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            aa[n * i + j] = a[i + 1][j + 1];
        }
    }
    for (i = 0; i < n + 1; i++) {
        delete[] a[i];
    }
    delete[] a;
    for (i = 0; i < n; i++) {
        b[i] = b[i + 1];
        x[i] = x[i + 1];
    }

    bool ppp;
    double nevv;
    ppp = SolveSystem(n, aa, b, x, kk, nevv);
    if (static_cast<int>(ppp) == 0) {
        cout << "determinant=0";
        delete[] aa;
        delete[] copb;
        delete[] b;
        delete[] x;
        fin.close();
        return -4;
    }
    double error = err(b, n);
    cout << "Residual: " << nevv << "\n";
    cout << "Error: " << error << '\n';
    cout << "Solution: "
         << " " << '\n';
    for (i = 0; i < m; i++) {
        printf("%.3e ", b[i]);
    }
    cout << "\n";
    delete[] copb;
    delete[] aa;
    delete[] b;
    delete[] x;
    fin.close();
    return 0;
}

/*bool Solve(int n, double *m, double *b, double *x, int *ind, int thread_num,
           int p, double *time)
{
    double elmax;
    int indmax;
    int s;
    double buf;
    double norma;
    int i;
    double time1 = currentTimeNano();
    for (i = 0; i < n; i++) {
        ind[i] = i;
    }
    norma = norm(n, m);
    for (int k = 0; k < n; k++) {
        // synchronize(p);
        if (thread_num == 0) {

            indmax = k;
            elmax = fabs(m[k * n + k]);
            for (int j = k + 1; j < n; j++) {
                if (fabs(m[k * n + j]) > elmax) {
                    elmax = fabs(m[k * n + j]);
                    indmax = j;
                }
            }

            // cout<<"indmax "<<indmax<<endl;
            // cout<<endl;
            // cout<<" ";
            s = ind[k];
            ind[k] = ind[indmax];
            // ind[k]=indmax;
            ind[indmax] = s;

            for (int i = 0; i < n; i++) {
                buf = m[i * n + k];
                m[i * n + k] = m[i * n + indmax];
                m[i * n + indmax] = buf;
            }

            if (fabs(m[k * n + k]) <= (const10)) {
                // cout<<"singular matrix";
                // error[0]=1;
                //
                // synchronize(p);
                // return 0;
                return true;
            }

            // cout<<"del"<<m[k*n+k];
            buf = 1.0 / m[k * n + k];
            // cout<<buf<<" ";
            for (int j = k; j < n; j++) {
                m[k * n + j] *= buf;
                // m[k*n+j]=m[k*n+j]/m[k*n+k];
                // cout<<"mainstr"<<m[k*n+j]<<" ";
            }
            // cout<<1.0 /m[k*n+k]<<'/n';
            // cout<<buf<<" ";
            b[k] *= buf;
            // cout<<b[k]<<" ";
            // cout<<endl;
        }
        synchronize(p);
        // cout<<b[k]<<'\n';
        for (int i = thread_num; i < k; i += p) {
            buf = m[i * n + k];
            for (int j = k; j < n; j++) {
                m[i * n + j] -= m[k * n + j] * buf;
            }
            b[i] -= b[k] * buf;
            // cout<<"thred no "<<thread_num<<" "<<buf<<endl;
        }
        // synchronize(p);
        for (int i = k + 1 + thread_num; i < n; i += p) {
            buf = m[i * n + k];
            for (int j = k; j < n; j++) {
                m[i * n + j] -= m[k * n + j] * buf;
                // m[i*n+j]=m[i*n+j]-buf*m[k*n+j];
                // cout<<"mem"<<m[i*n+j]<<" ";
            }
            // cout<<buf<<" ";
            b[i] -= b[k] * buf;
            // b[i]=b[i]-b[k]*buf;
            // cout<<"kek"<<b[i]<<" ";
            // cout<<"thred no "<<thread_num<<" "<<buf<<endl;
        }
        // synchronize(p);
        // cout<<"ans"<<b[k]<<" ";
    }
    // cout<<'\n';

    // for (int i=0;i<n;i++)
    // {
    //    cout<<ind[i]<<" ";
    // }
    // synchronize(p);
    // s++;
    if (thread_num == 0) {
        for (int i = 0; i < n; i++) {
            x[ind[i]] = b[i];
        }
        // synchronize(p);
    }

    // cout<<'\n';
    // cout<<x[3];
    // cout<<'\n';
    // cout<<norma;
    time1 = currentTimeNano() - time1;
    time[thread_num] = time1;

    return false;
}
*/
/*bool SolveSystem(int n, double *a, double *b, double *x, int *index,
                 int my_rank, int total_threads, double *time)
{
    int i;
    int j;
    int k;
    int first_row;
    int last_row;
    double tmp;
    double time1 = currentTimeNano();

    for (i = 0; i < n; i++) {
        index[i] = i;
    }

    for (i = 0; i < n; i++) {
        if (my_rank == 0) {
            k = i;
            for (j = i + 1; j < n; j++)
                if (fabs(a[i * n + k]) < fabs(a[i * n + j]))
                    k = j;

            j = index[i];
            index[i] = index[k];
            index[k] = j;

            for (j = 0; j < n; j++) {
                tmp = a[j * n + i];
                a[j * n + i] = a[j * n + k];
                a[j * n + k] = tmp;
            }

            tmp = a[i * n + i];
            if (abs(tmp) < eps) {
                return true;
            }
            tmp = 1.0 / tmp;
            for (j = i; j < n; j++)
                a[i * n + j] *= tmp;
            b[i] *= tmp;
        }
        synchronize(total_threads);

        first_row = (n - i - 1) * my_rank;
        first_row = first_row / total_threads + i + 1;
        last_row = (n - i - 1) * (my_rank + 1);
        last_row = last_row / total_threads + i + 1;

        for (j = first_row; j < last_row; j++) {
            tmp = a[j * n + i];
            for (k = i; k < n; k++)
                a[j * n + k] -= tmp * a[i * n + k];
            b[j] -= tmp * b[i];
        }
        synchronize(total_threads);
    }

    if (my_rank == 0)
        for (i = n - 1; i >= 0; i--) {
            tmp = b[i];
            for (j = i + 1; j < n; j++)
                tmp -= a[i * n + j] * x[index[j]];
            x[index[i]] = tmp;
        }
    time1 = currentTimeNano() - time1;
    time[my_rank] = time1;

    return false;
}
*/
/*
bool SolveSystem(int n, double *a, double *b, double *x, int my_rank,
                 int Threadcount, double *time)
{
    int i;
    int j;
    int k;
    int ii;
    int first_row;
    int last_row;
    double tmp;
    double time1 = currentTimeNano();
    for (i = 0; i < n; i++) {
        if (my_rank == 0) {
            for (ii = i; ii < n; ii++) {
                if (abs(a[ii * n + i]) > abs(a[i * n + i])) { //замена
                    for (k = i; k < n; k++) {

                        tmp = a[ii * n + k];
                        a[ii * n + k] = a[i * n + k];
                        a[i * n + k] = tmp;
                    }

                    tmp = b[ii];
                    b[ii] = b[i];
                    b[i] = tmp;
                }
            }
            tmp = a[i * n + i];
            if (abs(tmp) < eps) {
                return true;
            }
            tmp = 1.0 / tmp;
            for (j = i; j < n; j++) {
                a[i * n + j] *= tmp;
            }
            b[i] *= tmp;
        }
      //  synchronize(Threadcount);

        first_row = (n - i - 1) * my_rank;
        first_row = first_row / Threadcount + i + 1;
        last_row = (n - i - 1) * (my_rank + 1);
        last_row = last_row / Threadcount + i + 1;

        for (j = first_row; j < last_row; j++) {
            tmp = a[j * n + i];
            for (k = i; k < n; k++) {
                a[j * n + k] -= tmp * a[i * n + k];
            }
            b[j] -= tmp * b[i];
        }
        synchronize(Threadcount);
    }

    if (my_rank == 0) {
        for (i = n - 1; i >= 0; i--) {
            tmp = b[i];
            for (j = i + 1; j < n; j++) {
                tmp -= a[i * n + j] * x[j];
            }
            x[i] = tmp;
        }
    }
    time1 = currentTimeNano() - time1;
    time[my_rank] = time1;
    return false;
}
*/
/*bool SolveSystem(int n, double *a, double *b, double *x, int my_rank,
                 int Threadcount, double *time)
{
    int i;
    int j;
    int k;
    int ii;
    int first_row;
    int last_row;
    double tmp;
    double time1 = currentTimeNano();
    for (i = 0; i < n; i++) {
        if (my_rank == 0) {
            for (ii = i; ii < n; ii++) {
                if (abs(a[ii * n + i]) > abs(a[i * n + i])) { //замена
                    for (k = i; k < n; k++) {

                        tmp = a[ii * n + k];
                        a[ii * n + k] = a[i * n + k];
                        a[i * n + k] = tmp;
                    }

                    tmp = b[ii];
                    b[ii] = b[i];
                    b[i] = tmp;
                }
            }
            tmp = a[i * n + i];
            if (abs(tmp) < eps) {
                return true;
            }
            tmp = 1.0 / tmp;
            for (j = i; j < n; j++) {
                a[i * n + j] *= tmp;
            }
            b[i] *= tmp;
        }
      //  synchronize(Threadcount);

        first_row = (n - i - 1) * my_rank;
        first_row = first_row / Threadcount + i + 1;
        last_row = (n - i - 1) * (my_rank + 1);
        last_row = last_row / Threadcount + i + 1;

        for (j = first_row; j < last_row; j++) {
            tmp = a[j * n + i];
            for (k = i; k < n; k++) {
                a[j * n + k] -= tmp * a[i * n + k];
            }
            b[j] -= tmp * b[i];
        }
        synchronize(Threadcount);
    }

    if (my_rank == 0) {
        for (i = n - 1; i >= 0; i--) {
            tmp = b[i];
            for (j = i + 1; j < n; j++) {
                tmp -= a[i * n + j] * x[j];
            }
            x[i] = tmp;
        }
    }
    time1 = currentTimeNano() - time1;
    time[my_rank] = time1;
    return false;
}
*/