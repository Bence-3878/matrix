//
// Created by bence on 2023.05.20..
//

#ifndef HAZI_MATRIX_VEKTOR_H
#define HAZI_MATRIX_VEKTOR_H

#include <exception>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <sstream>

// Kivétel osztály az érvénytelen indexek kezeléséhez
class InvalidIndexException : public std::exception {
private:
    int invalidIndex; // Az érvénytelen indexet tároló változó.
    std::string message; // A hibaüzenetet tároló változó.

public:
    // A kivétel üzenetének és az érvénytelen indexnek a megadása a konstruktorban.
    InvalidIndexException(const char* message, int invalidIndex)
            : message(message), invalidIndex(invalidIndex) {}

    // Visszaadja az érvénytelen indexet.
    int getInvalidIndex() const {
        return invalidIndex;
    }

    // A hibaüzenet lekérdezése.
    const char* what() const noexcept override {
        return message.c_str();
    }
};

// Kivétel osztály a dimenziók nem megfelelősségének kezeléséhez
class DimensionMismatchException : public std::logic_error {
private:
    int dimension1;
    int dimension2;

public:
    DimensionMismatchException(const std::string& message, int dim1, int dim2)
            : std::logic_error(message), dimension1(dim1), dimension2(dim2) {}

    int getDimension1() const { return dimension1; }
    int getDimension2() const { return dimension2; }

    const char* what() const noexcept override {
        std::ostringstream oss;
        oss << std::logic_error::what()
            << " Első dimenzió: " << dimension1
            << ", második dimenzió: " << dimension2;
        return oss.str().c_str();
    }
};

template <typename T = int>
class Vektor {
    unsigned dimension;
    T* block;

public:
    // Alapértelmezett konstruktor
    Vektor();

    // Konstruktor a dimenzióval és alapértelmezett értékkel
    Vektor(unsigned dimension, const T& defaultValue);

    // Konstruktor a dimenzióval és nullákkal
    Vektor(unsigned dimension);

    // Másoló konstruktor
    Vektor(const Vektor<T>& vektor);

    // Konstruktor két vektor összefűzésével
    Vektor(const Vektor<T>& a, const Vektor<T>& b);

    // Destruktor
    ~Vektor();

    // Getter függvény a dimenzió lekérdezéséhez
    inline unsigned getDimension() const;

    // Elemek cseréje az adott indexeken
    void swap(int i, int j);

    // Az adott tartományból készít egy részvektort
    Vektor<T> getSubVector(int i, int j);

    // Ellenőrzi, hogy a vektor minden eleme nulla-e
    bool isVectorZero();

    // Getter függvény a blokk adatok lekérdezéséhez
    inline T* getData() const { return block; }

    // Dinamikus átméretezés
    void resize(int a);

    // Összeadás két vektor között
    Vektor<T> operator+(const Vektor<T>& vektor) const;

    // Kivonás két vektor között
    Vektor<T> operator-(const Vektor<T>& vektor) const;

    // Vektorok skaláris szorzata
    T operator*(const Vektor<T>& vektor) const;

    // Vektor skalárral való szorzása
    Vektor<T> operator*(const T& scalar) const;

    // Vektor skalárral való osztása
    Vektor<T> operator/(const T& scalar) const;

    // Indexelés (const változat)
    T& operator[](int i);

    // Indexelés (const változat)
    T& operator[](int i) const;

    // Értékadás operátor túlterhelése
    Vektor<T>& operator=(const Vektor<T>& vektor);

    // Globális << operátor túlterhelése (kiírás)
    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const Vektor<U>&);

    // Fájlba írás
    void writeToFile(const std::string& filename);

    // Fájlból olvasás
    void readFromFile(const std::string& filename);
};



// Alapértelmezett konstruktor
template<typename T>
Vektor<T>::Vektor() : dimension(0), block(nullptr) { }

// Konstruktor a dimenzióval és alapértelmezett értékkel
template<typename T>
Vektor<T>::Vektor(unsigned dimension, const T& defaultValue)
        : dimension(dimension) {
    block = new T[dimension];
    for (int i = 0; i < dimension; ++i)
        block[i] = defaultValue;
}

// Konstruktor a dimenzióval és nullákkal
template<typename T>
Vektor<T>::Vektor(unsigned dimension)
        : dimension(dimension) {
    block = new T[dimension];
    for (int i = 0; i < dimension; ++i)
        block[i] = 0;
}

// Másoló konstruktor
template<typename T>
Vektor<T>::Vektor(const Vektor<T>& vektor)
        : dimension(vektor.getDimension()) {
    block = new T[dimension];
    for (int i = 0; i < dimension; ++i)
        block[i] = vektor[i];
}

// Konstruktor két vektor összefűzésével
template<typename T>
Vektor<T>::Vektor(const Vektor<T>& a, const Vektor<T>& b)
        : dimension(a.getDimension() + b.getDimension()) {
    block = new T[dimension];
    int i;
    for (i = 0; i < a.getDimension(); ++i)
        block[i] = a[i];
    for (; i < getDimension(); ++i) {
        block[i] = b[i - a.getDimension()];
    }
}

// Destruktor
template<typename T>
Vektor<T>::~Vektor() {
    delete[] block;
}

// Getter függvény a dimenzió lekérdezéséhez
template<typename T>
inline unsigned Vektor<T>::getDimension() const {
    return dimension;
}

// Elemek cseréje az adott indexeken
template<typename T>
void Vektor<T>::swap(int i, int j) {
    if (i == j)
        return;
    T s = block[i];
    block[i] = block[j];
    block[j] = s;
}

// Az adott tartományból készít egy részvektort
template<typename T>
Vektor<T> Vektor<T>::getSubVector(int i, int j) {
    if (i > j)
        throw InvalidIndexException(
                "Kezdő cím nagyobb mint a vég cím", i);
    if (i < 0)
        throw InvalidIndexException(
                "Túlindexelés: az index nem lehet negatív", i);
    if (j >= dimension)
        throw InvalidIndexException(
                "Túlindexelés: az index meghaladja a dimenziót", j);

    Vektor<T> vektor(j - i + 1);
    for (int k = 0; k < j - i + 1; ++k)
        vektor[k] = block[i + k];

    return vektor;
}

// Ellenőrzi, hogy a vektor minden eleme nulla-e
template<typename T>
bool Vektor<T>::isVectorZero() {
    for (int i = 0; i < dimension; ++i)
        if (block[i] != 0)
            return false;
    return true;
}

template<typename T>
void Vektor<T>::resize(int a) {
    T* temp = new T[a];  // Új tömb létrehozása az új mérettel
    for (int i = 0; i < a && i < dimension; ++i)
        temp[i] = block[i];  // A meglévő elemek másolása az új tömbbe
    for (int i = dimension; i < a; ++i)
        temp[i] = 0;  // Az új elemek inicializálása 0-val
    dimension = a;  // Az új méret beállítása
    delete[] block;  // Az eredeti tömb felszabadítása
    block = temp;  // Az új tömbre mutató pointer beállítása
}

// Indexelés (nem const változat)
template<typename T>
T& Vektor<T>::operator[](int i) {
    if (i >= dimension)
        throw InvalidIndexException(
                "Túlindexelés: az index meghaladja a dimenziót", i);
    return block[i];
}

// Indexelés (const változat)
template<typename T>
T& Vektor<T>::operator[](int i) const {
    if (i >= dimension)
        throw InvalidIndexException(
                "Túlindexelés: az index meghaladja a dimenziót", i);
    return block[i];
}

// Értékadás operátor túlterhelése
template<typename T>
Vektor<T>& Vektor<T>::operator=(const Vektor<T>& vektor) {
    if (this != &vektor) {
        delete[] block;
        dimension = vektor.dimension;
        block = new T[dimension];
        for (unsigned i = 0; i < dimension; ++i)
            block[i] = vektor.block[i];
    }
    return *this;
}

// Globális szorzás operátor túlterhelése (T * Vektor)
template<typename T>
inline Vektor<T> operator*(const T& a, const Vektor<T>& vektor) {
    return vektor * a;  // Felhasználja a Vektor osztály belső szorzás operátorát
}

// Szorzás operátor túlterhelése (Vektor * T)
template<typename T>
Vektor<T> Vektor<T>::operator*(const T& a) const {
    Vektor<T> vektor(getDimension());
    for (int i = 0; i < getDimension(); ++i)
        vektor[i] = block[i] * a;
    return vektor;
}

// Osztás operátor túlterhelése (Vektor / T)
template<typename T>
Vektor<T> Vektor<T>::operator/(const T& a) const {
    Vektor<T> vektor(getDimension());
    for (int i = 0; i < getDimension(); ++i)
        vektor[i] = block[i] / a;
    return vektor;
}

// Szorzás operátor túlterhelése (Vektor * Vektor)
template<typename T>
T Vektor<T>::operator*(const Vektor<T>& vektor) const {
    if (getDimension() != vektor.getDimension())
        throw DimensionMismatchException(
                "Különböző vektor dimenziók: a művelet nem végrehajtható",
                getDimension(), vektor.getDimension());
    T szum = 0;
    for (int i = 0; i < getDimension(); ++i)
        szum = szum + block[i] * vektor[i];
    return szum;
}

// Kivonás operátor túlterhelése (Vektor - Vektor)
template<typename T>
Vektor<T> Vektor<T>::operator-(const Vektor<T>& vektor) const {
    if (getDimension() != vektor.getDimension())
        throw DimensionMismatchException(
                "Különböző vektor dimenziók: a művelet nem végrehajtható",
                getDimension(), vektor.getDimension());
    Vektor<T> vektor1(getDimension());
    for (int i = 0; i < getDimension(); ++i)
        vektor1[i] = block[i] - vektor[i];
    return vektor1;
}

// Összeadás operátor túlterhelése (Vektor + Vektor)
template<typename T>
Vektor<T> Vektor<T>::operator+(const Vektor<T>& vektor) const {
    if (getDimension() != vektor.getDimension())
        throw DimensionMismatchException(
                "Különböző vektor dimenziók: a művelet nem végrehajtható",
                getDimension(), vektor.getDimension());
    Vektor<T> vektor1(getDimension());
    for (int i = 0; i < getDimension(); ++i)
        vektor1[i] = block[i] + vektor[i];
    return vektor1;
}

// Globális << operátor túlterhelése (kiírás)
template<typename U>
std::ostream& operator<<(std::ostream& ostream, const Vektor<U>& vektor) {
    for (int i = 0; i < vektor.getDimension() - 1; ++i)
        ostream << vektor[i] << ", ";
    ostream << vektor[vektor.getDimension() - 1];

    return ostream;
}

// Fájlba írás
template<typename T>
void Vektor<T>::writeToFile(const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Hiba a fájl megnyitásakor: " + filename);
    }

    try {
        // Dimenzió írása a fájlba
        file.write(reinterpret_cast<const char*>(&dimension), sizeof(unsigned));

        // Vektor elemeinek írása a fájlba
        file.write(reinterpret_cast<const char*>(block), dimension * sizeof(T));

        file.close();
    } catch (const std::exception& e) {
        // Hibakezelés: Kivétel dobódott a fájlba történő írás során
        throw std::runtime_error("Hiba a fájlba történő írás során: " +
                                 std::string(e.what()));
    }
}

// Fájlból olvasás
template<typename T>
void Vektor<T>::readFromFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Hiba a fájl megnyitásakor: " + filename);
    }

    try {
        // Dimenzió beolvasása a fájlból
        file.read(reinterpret_cast<char*>(&dimension), sizeof(unsigned));

        // Dinamikus memóriaterület foglalása a vektor elemeinek
        block = new T[dimension];

        // Vektor elemeinek beolvasása a fájlból
        file.read(reinterpret_cast<char*>(block), dimension * sizeof(T));

        file.close();
    } catch (const std::exception& e) {
        // Hibakezelés: Kivétel dobódott a fájlból történő olvasás során
        delete[] block;  // Töröljük a már foglalt memóriaterületet
        block = nullptr; // Nullpointerre állítjuk a blokkot
        throw std::runtime_error("Hiba a fájlból történő olvasás során: " +
                                 std::string(e.what()));
    }
}

#endif // HAZI_MATRIX_VEKTOR_H
