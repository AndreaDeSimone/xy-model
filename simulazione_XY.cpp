#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <tuple>
#include <fstream>
#include <numeric>

#include <filesystem>
#include <cstdio>

// Ensure M_PI is available on MSVC if not provided by <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
namespace fs = std::filesystem;

// ================= PARAMETRI =================
const int L = 16;
double T = 1.5;
const double T_min = 0.1;
const double intervallo_T = 0.1;
// ================= DEFINIZIONI BASE =================
double beta_ = 1.0 / T;
const int N = L * L;
const int punti = (T - T_min) / intervallo_T;
// ================= RANDOM =================
std::mt19937 rng(std::random_device{}());
std::uniform_real_distribution<double> uni01(0.0, 1.0);
std::uniform_real_distribution<double> uni_angle(0.0, 2.0 * M_PI);
std::uniform_int_distribution<int> site_dist(0, L - 1);

using Matrix = vector<vector<double>>;

/* NEL CASO MALSANO CI VENISSE IN MENTE DI CAMBIARE LA DINAMICA
// ================= INTERAZIONE =================
vector<tuple<int, int, double>> Span;
const int raggio_pv = 1;

double metrica(int i, int j) {
    return sqrt(static_cast<double>(i * i + j * j));
}
// modulo sempre non-negativo (equivalente al comportamento di % in Python)
inline int pmod(int a, int m) {
    return ((a % m) + m) % m;
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

vector<tuple<int, int, double>>
double H_local(const Matrix& M, int i, int j) {
    double H = 0.0;
    for (auto& s : Span) {
        int di, dj; double dist;
        tie(di, dj, dist) = s;
        int ni = pmod(i + di, L);
        int nj = pmod(j + dj, L);
        H += -cos(M[i][j] - M[ni][nj]);
        }
    }
    return H;
}
*/

const int di[] = { 1, -1,  0,  0 };
const int dj[] = { 0,  0,  1, -1 };

// ================= DINAMICA =================
double H_local(const Matrix& M, int i, int j, double phi_2) {
    double H = 0.0;

    for (int m = 0; m < 4; m++)
    {
        int ni = (i + di[m] + L) % L;
        int nj = (j + dj[m] + L) % L;

        H += -cos(phi_2 - M[ni][nj]);
    }

    return H;
}

double H_tot(const Matrix& M) {
    double H = 0.0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            H += H_local(M, i, j, M[i][j]);
    return H / 4.0;
}

void dynamics(Matrix& M)
{
    // Estrai a caso un sito tra 0 e L-1
    int wi = site_dist(rng);
    int wj = site_dist(rng);

    // Nuovo angolo tra [0, 2pi)
    double phi_2 = uni_angle(rng);
    double phi_1 = M[wi][wj];
    double delta_E = H_local(M, wi, wj, phi_2) - H_local(M, wi, wj, phi_1);

    // Algoritmo di Metropolis (Bilancio dettagliato)
    if (delta_E <= 0.0 || uni01(rng) < exp(-beta_ * delta_E)) {
        M[wi][wj] = phi_2;
    }
}

// ================= MICROSTATO =================
double Ex_local(const Matrix& M, int i, int j) {
    double Ex = 0.0;
    Ex = cos(M[i][j] - M[i][(j + 1) % L]);
    return Ex;
}

double Ix_local(const Matrix& M, int i, int j) {
    double Ix = 0.0;
    Ix = sin(M[i][j] - M[i][(j + 1) % L]);
    return Ix;
}

double Ix_tot(const Matrix& M) {
    double Ix = 0.0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            Ix += Ix_local(M, i, j);
    return Ix;
}

double Ex_tot(const Matrix& M) {
    double Ex = 0.0;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            Ex += Ex_local(M, i, j);
    return Ex;
}

// ================= CHECKPOINT (salvataggio/caricamento matrice per T) =================
string cartellaMatrici() {
    return "matrici_L" + to_string(L);
}

// Costruisce il nome file associato a una temperatura, es. T=1.23 -> "matrici/matrice_T1.23.txt"
string matrixFilename(double T) {
    if (fabs(T) < 1e-9) T = 0.0; // evita il caso "-0.00"
    char buf[128];
    snprintf(buf, sizeof(buf), "%s/matrice_T%.2f.txt", cartellaMatrici().c_str(), T);
    return string(buf);
}

void saveMatrix(const Matrix& M, const string& filename) {
    ofstream out(filename);
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            out << M[i][j];
            if (j < L - 1) out << " ";
        }
        out << "\n";
    }
}

bool loadMatrix(Matrix& M, const string& filename) {
    ifstream in(filename);
    if (!in.is_open()) return false;
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            if (!(in >> M[i][j])) return false;
    return true;
}

// ================= SIMULAZIONE =================

template<typename... Args>
void print(const Args&... args) {
    (std::cout << ... << args) << '\n';
}

int main() {
    fs::create_directories(cartellaMatrici());  // crea "matrici_L128" se non esiste

    // generazione di una matrice casuale iniziale
    Matrix M(L, vector<double>(L));
    for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++)
            M[i][j] = 2 * M_PI * uni01(rng);

    vector<double> list_E, list_Ix, list_Ex, list_rho_s, list_T;

    // PUNTI DEL GRAFICO (temperatura decrescente)
    for (int i = 0; i < punti; i++) {
        string filename = matrixFilename(T);

        if (loadMatrix(M, filename)) {
            cout << "T=" << T << ": trovata matrice salvata, la carico e continuo la simulazione." << endl;
        }
        else {
            cout << "T=" << T << ": nessuna matrice salvata, parto da quella corrente." << endl;
            long long mc_1step = (long long)N * 1e4;
            for (long long step = 0; step < mc_1step; step++) {
                dynamics(M);
            }
        }
        double E_mean = 0.0;
        double Ix_mean = 0.0;
        double rho_s_mean = 0.0;
        int counter_mean = 0;
 

        T = T - intervallo_T;
        beta_ = 1.0 / T;
        cout << i << endl;

        int cooldown_ = 0;

        // STEP MONTECARLO
        long long mc_step = (long long)N * 1e4;
        for (long long step = 0; step < mc_step; step++) {
            dynamics(M);


            // MISURE
            if (step % N == 0) {
                double E = H_tot(M);
                double Ix = Ix_tot(M);
                double Ex = Ex_tot(M);
                double rho_s = (Ex - beta_ * (Ix * Ix)) / N;
                E_mean = E_mean + E;
                Ix_mean = Ix_mean + Ix;
                rho_s_mean = rho_s_mean + rho_s;
                counter_mean++;
            }
        }


        list_E.push_back(E_mean / (counter_mean * N));
        list_Ix.push_back(Ix_mean / (counter_mean * N));
        list_rho_s.push_back(rho_s_mean / counter_mean );
        print(counter_mean);
        print(H_tot(M)/N);
        print(Ix_tot(M)/N);
        print(rho_s_mean / counter_mean);
        counter_mean = 0;
        list_T.push_back(T);
        saveMatrix(M, filename);
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
