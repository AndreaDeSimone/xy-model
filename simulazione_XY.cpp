#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <tuple>
#include <fstream>
#include <numeric>

// Ensure M_PI is available on MSVC if not provided by <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

// ================= PARAMETRI =================
const int sqrtN = 16;
double T = 1.5;
const double T_min = 0.5;
const int punti = 50;
const bool Condizioni_al_bordo_periodiche = true;
// ================= DEFINIZIONI BASE =================
double beta_ = 1.0 / T;
const int N = sqrtN * sqrtN;
const double intervallo_T = (T - T_min) / (punti - 1);
// ================= RANDOM =================
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> uni01(0.0, 1.0);
std::normal_distribution<double> gauss(0.0, 1.0);
std::uniform_int_distribution<int> site_dist(0, sqrtN - 1);

using Matrix = vector<vector<double>>;

// ================= INTERAZIONE =================
const int raggio_pv = 1;

double metrica(int i, int j) {
    return sqrt(static_cast<double>(i * i + j * j));
}

// Trova le posizioni dei vicini entro una certa distanza (analogo di Span() in Python)
vector<tuple<int, int, double>> buildSpan(int raggio_pv) {
    vector<tuple<int, int, double>> lista;
    for (int i = -raggio_pv; i <= raggio_pv; i++) {
        for (int j = -raggio_pv; j <= raggio_pv; j++) {
            double d = metrica(i, j);
            if (d <= raggio_pv && i != j) {
                lista.push_back(make_tuple(i, j, d));
            }
        }
    }
    return lista;
}

vector<tuple<int, int, double>> Span;

// modulo sempre non-negativo (equivalente al comportamento di % in Python)
inline int pmod(int a, int m) {
    return ((a % m) + m) % m;
}

// ================= DINAMICA =================
double H_local(const Matrix& M, int i, int j) {
    double H = 0.0;
    for (auto& s : Span) {
        int di, dj; double dist;
        tie(di, dj, dist) = s;
        // NB: questa condizione è sempre vera perché Condizioni_al_bordo_periodiche == true
        // (esattamente come nel codice Python originale, dove "or True" annulla il controllo)
        if ((0 <= i + di && i + di < sqrtN && 0 <= j + dj && j + dj < sqrtN) ||
            Condizioni_al_bordo_periodiche) {
            int ni = pmod(i + di, sqrtN);
            int nj = pmod(j + dj, sqrtN);
            H += -cos(M[i][j] - M[ni][nj]);
        }
    }
    return H;
}

double Ex_local(const Matrix& M, int i, int j) {
    double Ex = 0.0;
    for (auto& s : Span) {
        int di, dj; double dist;
        tie(di, dj, dist) = s;
        if ((0 <= i + di && i + di < sqrtN && 0 <= j + dj && j + dj < sqrtN) ||
            Condizioni_al_bordo_periodiche) {
            if (di == 0) { // asse x
                int ni = pmod(i + di, sqrtN);
                int nj = pmod(j + dj, sqrtN);
                Ex += cos(M[i][j] - M[ni][nj]);
            }
        }
    }
    return Ex;
}

double Ix_local(const Matrix& M, int i, int j) {
    double Ix = 0.0;
    for (auto& s : Span) {
        int di, dj; double dist;
        tie(di, dj, dist) = s;
        if ((0 <= i + di && i + di < sqrtN && 0 <= j + dj && j + dj < sqrtN) ||
            Condizioni_al_bordo_periodiche) {
            if (di == 0) { // asse x
                int ni = pmod(i + di, sqrtN);
                int nj = pmod(j + dj, sqrtN);
                Ix += sin(M[i][j] - M[ni][nj]);
            }
        }
    }
    return Ix;
}

double H_tot(const Matrix& M) {
    double H = 0.0;
    for (int i = 0; i < sqrtN; i++)
        for (int j = 0; j < sqrtN; j++)
            H += H_local(M, i, j);
    return H / 4.0;
}

double Ix_tot(const Matrix& M) {
    double Ix = 0.0;
    for (int i = 0; i < sqrtN; i++)
        for (int j = 0; j < sqrtN; j++)
            Ix += Ix_local(M, i, j);
    return Ix / 2.0;
}

double Ex_tot(const Matrix& M) {
    double Ex = 0.0;
    for (int i = 0; i < sqrtN; i++)
        for (int j = 0; j < sqrtN; j++)
            Ex += Ex_local(M, i, j);
    return Ex / 2.0;
}

template<typename... Args>
void print(const Args&... args) {
    (std::cout << ... << args) << '\n';
}

int main() {
    Span = buildSpan(raggio_pv);

    // stampa lo Span, come il print(Span) del codice Python originale
    cout << "Span: ";
    for (auto& s : Span) {
        int di, dj; double dist;
        tie(di, dj, dist) = s;
        cout << "(" << di << ", " << dj << ", " << dist << ") ";
    }
    cout << endl;

    // generazione di una matrice casuale iniziale
    Matrix M(sqrtN, vector<double>(sqrtN));
    for (int i = 0; i < sqrtN; i++)
        for (int j = 0; j < sqrtN; j++)
            M[i][j] = 2 * M_PI * uni01(rng);

    const int passi = 50000;
    const int passo_iniziale = 20 * passi;
    const int cooldown = 10;

    vector<double> list_E, list_Ix, list_Ex, list_rho_s, list_T;

    // PUNTI DEL GRAFICO (temperatura decrescente)
    for (int i = 0; i < punti; i++) {
        double E = H_tot(M);
        double Ix = Ix_tot(M);
        double Ex = Ex_tot(M);

        double E_mean = 0.0;
        double Ix_mean = 0.0;
        double rho_s_mean = 0.0;
        int counter_mean = 0;

        T = T - intervallo_T;
        beta_ = 1.0 / T;
        cout << i << endl;

        int passo = (i == 0) ? passo_iniziale : passi;

        int cooldown_ = 0;

        // STEP MONTECARLO
        for (int step = 0; step < passo; step++) {
            int wi = site_dist(rng);
            int wj = site_dist(rng);

            // spostamento casuale dell'angolo nel sito scelto
            double x = sqrt(M_PI) * gauss(rng) / 2.0;

            Matrix M_new = M; // solo il sito (wi,wj) verrà modificato
            double val = fmod(M[wi][wj] + x, 2 * M_PI);
            if (val < 0) val += 2 * M_PI;
            M_new[wi][wj] = val;

            double delta_E = H_local(M_new, wi, wj) - H_local(M, wi, wj);

            // Bilancio dettagliato
            if (delta_E <= 0 || uni01(rng) < exp(-beta_ * delta_E)) {
                double old_Ix_local = Ix_local(M, wi, wj);
                double old_term_for_Ex = Ex_local(M, wi, wj);
                double new_Ix_local = Ix_local(M_new, wi, wj);
                double new_Ex_local = Ex_local(M_new, wi, wj);

                M = M_new;
                E = E + delta_E;
                Ix = Ix + new_Ix_local - old_Ix_local;
                Ex = Ex + new_Ex_local - old_term_for_Ex;
            }

            double rho_s = (Ex - beta_ * (Ix * Ix)) / N;

            // si iniziano a prendere le misure dopo il primo 90% dei passi
            if (passo - step <= 0.1 * passi) {
                if (cooldown_ >= cooldown) {
                    E_mean = E;
                    Ix_mean = Ix;
                    rho_s_mean = rho_s;
                    counter_mean++;
                    cooldown_ = 0;
                }
                else {
                    cooldown_++;
                }
            }
        }


        list_E.push_back(E_mean / (counter_mean * N));
        list_Ix.push_back(Ix_mean / (counter_mean * N));
        list_rho_s.push_back(rho_s_mean / counter_mean );
        print(counter_mean);
        counter_mean = 0;
        list_T.push_back(T);
    }

    // salvataggio risultati su CSV (al posto dei grafici)
    ofstream out("risultati_xy.csv");
    out << "T,E,Ix,rho_s\n";
    for (size_t i = 0; i < list_T.size(); i++) {
        out << list_T[i] << "," << list_E[i] << "," << list_Ix[i] << "," << list_rho_s[i] << "\n";
    }
    out.close();

    cout << "Simulazione completata. Risultati salvati in risultati_xy.csv" << endl;

    return 0;


}
