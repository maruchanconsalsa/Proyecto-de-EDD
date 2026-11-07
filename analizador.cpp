// analizador.cpp - Analizador de movimiento con polinomios
// Estructuras usadas (nodos y punteros):
//   - Lista enlazada ordenada: Polinomio
//   - Pila: historial para deshacer
//   - Cola: bitacora de operaciones
//
// Modelo fisico: s(t) posicion, v(t) = s'(t), a(t) = s''(t)
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstdlib>
#include <algorithm>
using namespace std;


// PARTE 1: PILA 
template <typename T>
class Pila {
private:
    struct NodoP {
        T dato;
        NodoP* sig;
        NodoP(const T& d, NodoP* s) : dato(d), sig(s) {}
    };
    NodoP* tope;
    int n;

    // Copia prohibida 
    Pila(const Pila&);
    Pila& operator=(const Pila&);

public:
    Pila() : tope(NULL), n(0) {}
    ~Pila() { T basura; while (pop(basura)) {} }

    void push(const T& dato) { tope = new NodoP(dato, tope); n++; }

    bool pop(T& salida) {
        if (!tope) return false;
        NodoP* aux = tope;
        salida = aux->dato;
        tope = aux->sig;
        delete aux;
        n--;
        return true;
    }

    bool vacia() const { return tope == NULL; }
    int tamano() const { return n; }
};

// PARTE 2: COLA
template <typename T>
class Cola {
private:
    struct NodoC {
        T dato;
        NodoC* sig;
        explicit NodoC(const T& d) : dato(d), sig(NULL) {}
    };
    NodoC* frente;
    NodoC* fin;
    int n;

    // Copia prohibida 
    Cola(const Cola&);
    Cola& operator=(const Cola&);

public:
    Cola() : frente(NULL), fin(NULL), n(0) {}
    ~Cola() { T basura; while (desencolar(basura)) {} }

    void encolar(const T& dato) {
        NodoC* nuevo = new NodoC(dato);
        if (fin) fin->sig = nuevo; else frente = nuevo;
        fin = nuevo;
        n++;
    }

    bool desencolar(T& salida) {
        if (!frente) return false;
        NodoC* aux = frente;
        salida = aux->dato;
        frente = frente->sig;
        if (!frente) fin = NULL;
        delete aux;
        n--;
        return true;
    }

    // Lee la posicion pos (0 = la mas antigua) sin sacarla
    bool verEn(int pos, T& salida) const {
        if (pos < 0 || pos >= n) return false;
        NodoC* p = frente;
        for (int i = 0; i < pos; i++) p = p->sig;
        salida = p->dato;
        return true;
    }

    bool vacia() const { return frente == NULL; }
    int tamano() const { return n; }
};

// PARTE 3: LISTA ENLAZADA ORDENADA - Polinomio
// Un termino: coef * t^exp
struct Nodo {
    double coef;
    int exp;
    Nodo* sig;
    Nodo(double c, int e) : coef(c), exp(e), sig(NULL) {}
};

class Polinomio {
private:
    Nodo* cabeza;                       // termino de mayor exponente
    void copiarDe(const Polinomio& o);

public:
    Polinomio();
    Polinomio(const Polinomio& o);
    Polinomio& operator=(const Polinomio& o);
    ~Polinomio();

    // Devuelve true si combino con un termino semejante
    bool insertarTermino(double coef, int exp);
    void limpiar();
    bool estaVacio() const { return cabeza == NULL; }
    int grado() const;                  // -1 si es el polinomio cero
    int cantidadTerminos() const;

    Polinomio sumar(const Polinomio& o) const;
    Polinomio restar(const Polinomio& o) const;
    Polinomio multiplicar(const Polinomio& o) const;

    Polinomio derivar() const;
    Polinomio integrar() const;         // con C = 0
    double evaluar(double t) const;     // Horner para terminos dispersos

    std::vector<double> raices(double a, double b) const;
    std::string aTexto(char var = 't') const;
};

// PARTE 3 : implementacion de Polinomio
namespace {
    const double EPS = 1e-12;           // tolerancia para considerar cero

    // Evita guardar la misma raiz dos veces seguidas
    void agregarRaiz(vector<double>& res, double x) {
        if (res.empty() || fabs(res.back() - x) > 1e-6) res.push_back(x);
    }
}

// ---------------------------- Memoria ----------------------------
Polinomio::Polinomio() : cabeza(NULL) {}

Polinomio::Polinomio(const Polinomio& o) : cabeza(NULL) { copiarDe(o); }

Polinomio& Polinomio::operator=(const Polinomio& o) {
    if (this != &o) {
        limpiar();
        copiarDe(o);
    }
    return *this;
}

Polinomio::~Polinomio() { limpiar(); }

void Polinomio::limpiar() {
    while (cabeza) {
        Nodo* aux = cabeza;
        cabeza = cabeza->sig;
        delete aux;
    }
}

void Polinomio::copiarDe(const Polinomio& o) {
    Nodo* ultimo = NULL;
    for (Nodo* p = o.cabeza; p; p = p->sig) {
        Nodo* nuevo = new Nodo(p->coef, p->exp);
        if (!ultimo) cabeza = nuevo; else ultimo->sig = nuevo;
        ultimo = nuevo;
    }
}

// ---------------------------- Insercion ----------------------------
bool Polinomio::insertarTermino(double coef, int exp) {
    if (exp < 0) return false;

    // Buscar la posicion (orden descendente)
    Nodo* prev = NULL;
    Nodo* act = cabeza;
    while (act && act->exp > exp) { prev = act; act = act->sig; }

    // Mismo exponente: sumar coeficientes
    if (act && act->exp == exp) {
        act->coef += coef;
        if (fabs(act->coef) < EPS) {    // quedo en cero: se elimina
            if (prev) prev->sig = act->sig; else cabeza = act->sig;
            delete act;
        }
        return true;
    }

    // Exponente nuevo: crear nodo (no se guardan ceros)
    if (fabs(coef) < EPS) return false;
    Nodo* nuevo = new Nodo(coef, exp);
    nuevo->sig = act;
    if (prev) prev->sig = nuevo; else cabeza = nuevo;
    return false;
}

int Polinomio::grado() const { return cabeza ? cabeza->exp : -1; }

int Polinomio::cantidadTerminos() const {
    int c = 0;
    for (Nodo* p = cabeza; p; p = p->sig) c++;
    return c;
}

// ---------------------------- Algebra ----------------------------
Polinomio Polinomio::sumar(const Polinomio& o) const {
    Polinomio r(*this);
    for (Nodo* p = o.cabeza; p; p = p->sig) r.insertarTermino(p->coef, p->exp);
    return r;
}

Polinomio Polinomio::restar(const Polinomio& o) const {
    Polinomio r(*this);
    for (Nodo* p = o.cabeza; p; p = p->sig) r.insertarTermino(-p->coef, p->exp);
    return r;
}

Polinomio Polinomio::multiplicar(const Polinomio& o) const {
    Polinomio r;
    for (Nodo* a = cabeza; a; a = a->sig)
        for (Nodo* b = o.cabeza; b; b = b->sig)
            r.insertarTermino(a->coef * b->coef, a->exp + b->exp);
    return r;
}

// ---------------------------- Calculo ----------------------------
// d/dt (c t^e) = c*e t^(e-1)
Polinomio Polinomio::derivar() const {
    Polinomio r;
    for (Nodo* p = cabeza; p; p = p->sig)
        if (p->exp > 0) r.insertarTermino(p->coef * p->exp, p->exp - 1);
    return r;
}

// integral de c t^e = c/(e+1) t^(e+1)
Polinomio Polinomio::integrar() const {
    Polinomio r;
    for (Nodo* p = cabeza; p; p = p->sig)
        r.insertarTermino(p->coef / (p->exp + 1), p->exp + 1);
    return r;
}

// Horner: entre terminos se multiplica por t^(diferencia de exponentes)
double Polinomio::evaluar(double t) const {
    if (!cabeza) return 0.0;
    double res = cabeza->coef;
    int expPrev = cabeza->exp;
    for (Nodo* p = cabeza->sig; p; p = p->sig) {
        res = res * pow(t, (double)(expPrev - p->exp)) + p->coef;
        expPrev = p->exp;
    }
    return res * pow(t, (double)expPrev);
}

// Barrido en N tramos y biseccion donde cambia el signo
vector<double> Polinomio::raices(double a, double b) const {
    vector<double> res;
    if (!cabeza || grado() == 0 || a >= b) return res;

    const int N = 2000;
    const double h = (b - a) / N;
    double x0 = a, f0 = evaluar(a);
    if (fabs(f0) < EPS) agregarRaiz(res, x0);

    for (int i = 1; i <= N; i++) {
        double x1 = a + i * h, f1 = evaluar(x1);
        if (fabs(f1) < EPS) {
            agregarRaiz(res, x1);
        } else if (fabs(f0) >= EPS && f0 * f1 < 0) {
            double lo = x0, hi = x1, flo = f0;
            for (int k = 0; k < 100 && (hi - lo) > 1e-10; k++) {
                double mid = (lo + hi) / 2, fm = evaluar(mid);
                if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; }
            }
            agregarRaiz(res, (lo + hi) / 2);
        }
        x0 = x1; f0 = f1;
    }
    return res;
}

// ---------------------------- Texto ----------------------------
string Polinomio::aTexto(char var) const {
    if (!cabeza) return "0";
    ostringstream os;
    bool primero = true;
    for (Nodo* p = cabeza; p; p = p->sig) {
        double c = p->coef, ac = fabs(c);
        if (primero) { if (c < 0) os << "-"; }
        else os << (c < 0 ? " - " : " + ");

        // El coeficiente 1 se omite si hay variable
        if (fabs(ac - 1.0) > 1e-12 || p->exp == 0) os << ac;
        if (p->exp >= 1) os << var;
        if (p->exp > 1) os << "^" << p->exp;
        primero = false;
    }
    return os.str();
}

// PARTE 4: PROGRAMA PRINCIPAL
const int MAX_TERMINOS = 20;
const int MAX_EXP = 50;
const double MAX_COEF = 1e9;
const int MAX_BITACORA = 20;

// Lo que se guarda en la pila para deshacer
struct Registro {
    Polinomio copia;
    string descripcion;
};

// ---------------------- Lectura con validacion ----------------------
static string recortar(const string& s) {
    size_t i = s.find_first_not_of(" \t\r\n");
    if (i == string::npos) return "";
    size_t j = s.find_last_not_of(" \t\r\n");
    return s.substr(i, j - i + 1);
}

// Reemplaza std::to_string (no existe en Dev-C++)
static string aTextoEntero(int x) {
    ostringstream os;
    os << x;
    return os.str();
}

static string leerLinea(const string& msg) {
    cout << msg;
    string s;
    if (!getline(cin, s)) { cout << "\n[Fin de la entrada]\n"; exit(0); }
    return s;
}

// Acepta coma o punto decimal; rechaza texto, nan e inf
static bool aDouble(string s, double& out) {
    s = recortar(s);
    if (s.empty()) return false;
    replace(s.begin(), s.end(), ',', '.');
    if (s.find_first_not_of("0123456789+-.eE") != string::npos) return false;
    char* fin = NULL;
    double v = strtod(s.c_str(), &fin);
    if (fin == s.c_str() || *fin != '\0') return false;
    if (!(v == v) || fabs(v) > DBL_MAX) return false;
    out = v;
    return true;
}

static bool aEntero(const string& txt, int& out) {
    string s = recortar(txt);
    if (s.empty()) return false;
    char* fin = NULL;
    long v = strtol(s.c_str(), &fin, 10);
    if (fin == s.c_str() || *fin != '\0') return false;
    if (v < INT_MIN || v > INT_MAX) return false;
    out = (int)v;
    return true;
}

static int leerEntero(const string& msg, int lo, int hi) {
    while (true) {
        int v;
        if (aEntero(leerLinea(msg), v) && v >= lo && v <= hi) return v;
        cout << "  [Error] Ingrese un entero entre " << lo << " y " << hi << ".\n";
    }
}

static double leerDouble(const string& msg, double lo, double hi) {
    while (true) {
        double v;
        if (aDouble(leerLinea(msg), v) && v >= lo && v <= hi) return v;
        cout << "  [Error] Ingrese un numero real entre " << lo << " y " << hi << ".\n";
    }
}

static bool leerSiNo(const string& msg) {
    while (true) {
        string s = recortar(leerLinea(msg));
        if (s == "s" || s == "S") return true;
        if (s == "n" || s == "N") return false;
        cout << "  [Error] Responda con 's' o 'n'.\n";
    }
}

// ---------------------- Funciones del programa ----------------------
static string num(double x) {
    ostringstream os; os << setprecision(6) << x; return os.str();
}

// Anota en la cola; si pasa de 20 saca la mas antigua
static void registrar(Cola<string>& bitacora, const string& texto) {
    bitacora.encolar(texto);
    string descartado;
    if (bitacora.tamano() > MAX_BITACORA) bitacora.desencolar(descartado);
}

static Polinomio ingresarPolinomio(const string& nombre) {
    Polinomio p;
    cout << "\n--- Ingreso de " << nombre << " (forma: c*t^e) ---\n";
    int n = leerEntero("Cantidad de terminos (1-" + aTextoEntero(MAX_TERMINOS) + "): ",
                       1, MAX_TERMINOS);
    for (int i = 1; i <= n; i++) {
        cout << "Termino " << i << ":\n";
        double c = leerDouble("  Coeficiente: ", -MAX_COEF, MAX_COEF);
        int e = leerEntero("  Exponente (0-" + aTextoEntero(MAX_EXP) + "): ", 0, MAX_EXP);
        if (c == 0) {
            cout << "  [Aviso] Coeficiente 0: termino ignorado.\n";
            continue;
        }
        if (p.insertarTermino(c, e))
            cout << "  [Aviso] Exponente repetido: se combinaron terminos semejantes.\n";
    }
    cout << nombre << " = " << p.aTexto() << "\n";
    return p;
}

// Muestra el resultado y pregunta si se guarda en A (apila el estado anterior)
static void ofrecerGuardar(const Polinomio& r, const string& desc, Polinomio& A,
                           Pila<Registro>& historial, Cola<string>& bitacora) {
    cout << "Resultado: " << r.aTexto() << "\n";
    registrar(bitacora, desc + " => " + r.aTexto());
    if (leerSiNo("Guardar el resultado en A? (s/n): ")) {
        Registro reg; reg.copia = A; reg.descripcion = desc;
        historial.push(reg);
        A = r;
        cout << "A = " << A.aTexto() << "\n";
    }
}

static void buscarRaices(const Polinomio& p, const string& nombre, bool esVelocidad) {
    cout << "\n" << nombre << " = " << p.aTexto() << "\n";
    if (p.estaVacio()) { cout << "Es el polinomio cero: todo t es solucion.\n"; return; }
    if (p.grado() == 0) { cout << "Es una constante no nula: no tiene raices.\n"; return; }

    double a = leerDouble("Tiempo inicial (s, >= 0): ", 0, 1e6);
    double b = leerDouble("Tiempo final (s, > inicial): ", a + 1e-6, 1e6);

    vector<double> r = p.raices(a, b);
    if (r.empty()) {
        cout << "No se encontraron raices en [" << a << ", " << b << "].\n";
        return;
    }
    cout << fixed << setprecision(6);
    for (size_t i = 0; i < r.size(); i++)
        cout << "  t" << i + 1 << " = " << r[i] << " s\n";
    cout.unsetf(ios::fixed);

    if (esVelocidad) cout << "En esos instantes el movil se detiene o cambia de sentido.\n";
    else cout << "En esos instantes el movil pasa por el origen (s = 0).\n";
    cout << "Nota: raices dobles (sin cambio de signo) pueden no detectarse.\n";
}

static void mostrarMenu() {
    cout << "\n========= ANALIZADOR DE MOVIMIENTO CON POLINOMIOS =========\n"
         << "  1. Ingresar polinomio A (posicion s(t))\n"
         << "  2. Ingresar polinomio B\n"
         << "  3. Mostrar A y B\n"
         << "  4. A + B\n"
         << "  5. A - B\n"
         << "  6. A * B\n"
         << "  7. Derivar A (velocidad y aceleracion)\n"
         << "  8. Integrar A\n"
         << "  9. Evaluar A en un tiempo t (s, v, a)\n"
         << " 10. Raices de s(t) = 0\n"
         << " 11. Instantes de reposo v(t) = 0\n"
         << " 12. Deshacer ultimo cambio en A\n"
         << " 13. Ver bitacora de operaciones\n"
         << " 14. Cargar ejemplo en A: 4t^3 - 2t^2 + 5t - 7\n"
         << " 15. Salir\n";
}

int main() {
    Polinomio A, B;
    Pila<Registro> historial;           // estados anteriores de A
    Cola<string> bitacora;              // ultimas operaciones
    int op;

    do {
        mostrarMenu();
        op = leerEntero("Opcion: ", 1, 15);

        switch (op) {
            case 1: {
                Polinomio nuevo = ingresarPolinomio("A");
                Registro reg; reg.copia = A; reg.descripcion = "Ingreso de A";
                historial.push(reg);
                A = nuevo;
                registrar(bitacora, "Ingreso A = " + A.aTexto());
                break;
            }
            case 2:
                B = ingresarPolinomio("B");
                registrar(bitacora, "Ingreso B = " + B.aTexto());
                break;
            case 3:
                cout << "\nA = " << A.aTexto() << "  (grado " << A.grado() << ")\n"
                     << "B = " << B.aTexto() << "  (grado " << B.grado() << ")\n";
                break;
            case 4: ofrecerGuardar(A.sumar(B), "A + B", A, historial, bitacora); break;
            case 5: ofrecerGuardar(A.restar(B), "A - B", A, historial, bitacora); break;
            case 6: ofrecerGuardar(A.multiplicar(B), "A * B", A, historial, bitacora); break;
            case 7: {
                Polinomio v = A.derivar(), a = v.derivar();
                cout << "\ns(t) = " << A.aTexto() << "   [m]\n"
                     << "v(t) = " << v.aTexto() << "   [m/s]\n"
                     << "a(t) = " << a.aTexto() << "   [m/s^2]\n";
                ofrecerGuardar(v, "Derivar A (v = s')", A, historial, bitacora);
                break;
            }
            case 8:
                cout << "\n(Constante de integracion C = 0)\n";
                ofrecerGuardar(A.integrar(), "Integrar A", A, historial, bitacora);
                break;
            case 9: {
                double t = leerDouble("Tiempo t (s, >= 0): ", 0, 1e6);
                Polinomio v = A.derivar(), a = v.derivar();
                cout << "\nEn t = " << t << " s:\n"
                     << "  Posicion     s = " << num(A.evaluar(t)) << " m\n"
                     << "  Velocidad    v = " << num(v.evaluar(t)) << " m/s\n"
                     << "  Aceleracion  a = " << num(a.evaluar(t)) << " m/s^2\n";
                registrar(bitacora, "Evaluar A en t=" + num(t) + " => s=" + num(A.evaluar(t)));
                break;
            }
            case 10:
                buscarRaices(A, "s(t)", false);
                registrar(bitacora, "Busqueda de raices de s(t)");
                break;
            case 11:
                buscarRaices(A.derivar(), "v(t)", true);
                registrar(bitacora, "Busqueda de instantes de reposo v(t)=0");
                break;
            case 12: {
                Registro reg;
                if (historial.pop(reg)) {
                    A = reg.copia;
                    cout << "Se deshizo: " << reg.descripcion << ". A = " << A.aTexto() << "\n";
                    registrar(bitacora, "Deshacer: " + reg.descripcion);
                } else {
                    cout << "No hay cambios para deshacer.\n";
                }
                break;
            }
            case 13: {
                if (bitacora.vacia()) { cout << "La bitacora esta vacia.\n"; break; }
                cout << "\n--- Bitacora (de la mas antigua a la mas reciente) ---\n";
                string linea;
                for (int i = 0; i < bitacora.tamano(); i++) {
                    bitacora.verEn(i, linea);
                    cout << i + 1 << ". " << linea << "\n";
                }
                break;
            }
            case 14: {
                Registro reg; reg.copia = A; reg.descripcion = "Carga del ejemplo";
                historial.push(reg);
                A.limpiar();
                A.insertarTermino(4, 3); A.insertarTermino(-2, 2);
                A.insertarTermino(5, 1); A.insertarTermino(-7, 0);
                cout << "A = " << A.aTexto() << "\n";
                registrar(bitacora, "Cargar ejemplo A = " + A.aTexto());
                break;
            }
            case 15:
                cout << "Hasta luego.\n";
                break;
        }
    } while (op != 15);

    // Pausa para que la ventana de se cierre sola
    leerLinea("\nPresione Enter para cerrar...");
    return 0;
}
