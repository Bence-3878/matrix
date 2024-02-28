//
// Created by bence on 2023.05.20..
//

#ifndef HAZI_MATRIX_MATRIX_H
#define HAZI_MATRIX_MATRIX_H

#include "Vektor.h"
#include <fstream>
#include <exception>
#include <string>

// Kivételosztály a mátrixhoz
class MatrixException : public std::exception {
private:
    std::string message;

public:
    explicit MatrixException(const std::string& message) : message(message) {}

    virtual const char* what() const noexcept override {
        return message.c_str();
    }
};

// Típusparaméterként megadott elemekből (alapértelmezetten egész számokból)
// álló mátrixot definiáló osztály.
template <typename T = int>
class Matrix {
    // A mátrix sorainak és oszlopainak száma.
    unsigned row_count, column_count;
    // A mátrix elemeit tartalmazó blokk.
    Vektor<Vektor<T>> block;

public:
    // Alapértelmezett konstruktor, amely létrehoz egy üres mátrixot.
    Matrix();
    // Konstruktor, amely létrehoz egy adott méretű mátrixot.
    Matrix(int, int);
    // Konstruktor, amely létrehoz egy adott méretű mátrixot és feltölti a megadott tömbből.
    Matrix(int, int, T*, int);
    // Másoló konstruktor.
    Matrix(const Matrix<T>&);

    // Kiszámítja és visszaadja a mátrix determinánsát.
    T det() const;
    // Kiszámítja és visszaadja a mátrix inverzét.
    Matrix<T> inverse() const;
    // Visszaadja a mátrix egy részét.
    Matrix<T> getSubMatrix(int, int, int, int) const;
    // Kicseréli a mátrix két sorát.
    void swapRows(int, int);
    // Kicseréli a mátrix két oszlopát.
    void swapColumns(int, int);
    // Összefűzi a mátrixot egy másik mátrixszal vízszintesen.
    Matrix<T> mergeHorizontal(const Matrix<T>&) const;
    // Összefűzi a mátrixot egy másik mátrixszal függőlegesen.
    Matrix<T> mergeVertical(const Matrix<T>&) const;
    // Kiszámolja a mátrix transzponáltját.
    Matrix<T> transposition() const;
    // Elvégzi a mátrixon a Gauss-Jordan eliminációt.
    Matrix<T> GaussJordanElimination() const;
    // Kiírja a mátrix tartalmát egy fájlba.
    void writeToFile(const std::string& filename);
    // Beolvassa a mátrix tartalmát egy fájlból.
    void readFromFile(const std::string& filename);

    // Túlterhelt indexelő operátorok, amelyek lehetővé teszik a mátrix elemeinek elérését.
    Vektor<T> operator[](int) const;
    Vektor<T>& operator[](int);

    // Mátrixműveletek támogatása operátorok túlterhelésével.
    Matrix<T> operator-(const Matrix<T>&) const; // Mátrix kivonása
    Matrix<T> operator+(const Matrix<T>&) const; // Mátrix összeadása
    Matrix<T> operator*(const Matrix<T>&) const; // Mátrixok szorzása
    Matrix<T> operator*(const T&) const; // Skalárral való szorzás

    // Barát függvény, ami lehetővé teszi a mátrix kiírását egy kimeneti adatfolyamba.
    template<typename U>
    friend std::ostream& operator<<(std::ostream&, const Matrix<U>&);
};

// Alapértelmezett konstruktor, amely létrehoz egy üres mátrixot.
template<typename T>
Matrix<T>::Matrix() : row_count(0), column_count(0), block() {}

// Konstruktor, amely létrehoz egy adott méretű mátrixot.
template<typename T>
Matrix<T>::Matrix(int n, int m) : row_count(n), column_count(m), block(n, Vektor<T>(m)) {}

// Konstruktor, amely létrehoz egy adott méretű mátrixot és feltölti a megadott tömbből.
template<typename T>
Matrix<T>::Matrix(int n, int m, T* arr, int arr_size)
        : row_count(n), column_count(m), block(n, Vektor<T>(m)) {
    if (n * m < arr_size)
        throw std::domain_error("A tömb mérete nem egyezik a mátrix méretével.");

    int k = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (k < arr_size)
                block[i][j] = arr[k++];
            else
                block[i][j] = 0;
        }
    }
}

// Másoló konstruktor.
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix)
        : row_count(matrix.row_count), column_count(matrix.column_count),
          block(row_count, Vektor<T>(column_count)) {
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < column_count; ++j) {
            block[i][j] = matrix[i][j];
        }
    }
}

// Speciális függvény, amely egy egységmátrixot hoz létre.
template<typename T>
Matrix<T> unitMatrix(int size) {
    Matrix<T> matrix(size, size);
    for (int i = 0; i < size; ++i) {
        matrix[i][i] = 1;
    }
    return matrix;
}

// Kiszámítja és visszaadja a mátrix determinánsát.
template<typename T>
T Matrix<T>::det() const {
    if (column_count != row_count)
        throw std::domain_error("A mátrix nem négyzetes.");

    if (row_count == 2)
        return block[0][0] * block[1][1] - block[0][1] * block[1][0];

    try {
        T sum = block[0][0] * getSubMatrix(1, row_count - 1, 1, row_count - 1).det();
        for (int i = 1; i < row_count - 1; ++i) {
            if (i % 2 == 0)
                sum += block[i][0] * getSubMatrix(0, i - 1, 1, row_count - 1)
                        .mergeVertical(getSubMatrix(
                                i + 1, row_count - 1, 1, row_count - 1)).det();
            else
                sum -= block[i][0] * getSubMatrix(0, i - 1, 1, row_count - 1)
                        .mergeVertical(getSubMatrix(
                                i + 1, row_count - 1, 1, row_count - 1)).det();
        }
        if (row_count % 2 == 1)
            sum += block[row_count - 1][0] *
                    getSubMatrix(0, row_count - 2, 1, row_count - 1).det();
        else
            sum -= block[row_count - 1][0] *
                    getSubMatrix(0, row_count - 2, 1, row_count - 1).det();

        return sum;
    } catch (const std::exception& e) {
        throw std::runtime_error("Hiba a determináns számítása során: " + std::string(e.what()));
    }
}

// Kiszámítja és visszaadja a mátrix inverzét.
template<typename T>
Matrix<T> Matrix<T>::inverse() const {
    if (det() == 0)
        throw MatrixException("A mátrix nem invertálható, mert a determinánsa 0.");

    return mergeHorizontal(unitMatrix<T>(row_count))
            .GaussJordanElimination()
            .getSubMatrix(0, row_count - 1, row_count, 2 * row_count - 1);
}

// Visszaadja a mátrix egy részét.
template<typename T>
Matrix<T> Matrix<T>::getSubMatrix(int i, int j, int k, int l) const {
    if (i > j)
        throw std::out_of_range("A kezdő sorindex nagyobb, mint a vég sorindex.");
    if (i < 0 || i >= row_count)
        throw std::out_of_range("A kezdő sorindex nincs a megfelelő tartományban.");
    if (j >= row_count)
        throw std::out_of_range("A vég sorindex nincs a megfelelő tartományban.");
    if (k > l)
        throw std::out_of_range("A kezdő oszlopindex nagyobb, mint a vég oszlopindex.");
    if (k < 0 || k >= column_count)
        throw std::out_of_range("A kezdő oszlopindex nincs a megfelelő tartományban.");
    if (l >= column_count)
        throw std::out_of_range("A vég oszlopindex nincs a megfelelő tartományban.");

    Matrix<T> matrix(j - i + 1, l - k + 1);
    for (int i1 = 0; i1 < j - i + 1; ++i1) {
        matrix.block[i1] = block[i + i1].getSubVector(k, l);
    }
    return matrix;
}

// Kicseréli a mátrix két sorát.
template<typename T>
void Matrix<T>::swapRows(int i, int j) {
    block.swap(i, j);
}

// Kicseréli a mátrix két oszlopát.
template<typename T>
void Matrix<T>::swapColumns(int i, int j) {
    for (int k = 0; k < row_count; ++k) {
        block[k].swap(i, j);
    }
}

// Összefűzi a mátrixot egy másik mátrixszal vízszintesen.
template<typename T>
Matrix<T> Matrix<T>::mergeHorizontal(const Matrix<T>& matrix) const {
    if (row_count != matrix.row_count)
        throw MatrixException(
                "A mátrixok sorainak száma nem egyezik a horizontális összefűzésnél.");

    Matrix<T> result(row_count, column_count + matrix.column_count);
    for (int i = 0; i < row_count; ++i) {
        int j;
        for (j = 0; j < column_count; ++j) {
            result[i][j] = block[i][j];
        }
        for (; j < column_count + matrix.column_count; ++j) {
            result[i][j] = matrix[i][j - column_count];
        }
    }
    return result;
}

// Összefűzi a mátrixot egy másik mátrixszal függőlegesen.
template<typename T>
Matrix<T> Matrix<T>::mergeVertical(const Matrix<T>& matrix) const {
    if (column_count != matrix.column_count)
        throw MatrixException(
                "A mátrixok oszlopainak száma nem egyezik a függőleges összefűzésnél.");

    Matrix<T> result(row_count + matrix.row_count, column_count);
    for (int j = 0; j < column_count; ++j) {
        int i;
        for (i = 0; i < row_count; ++i) {
            result[i][j] = block[i][j];
        }
        for (; i < row_count + matrix.row_count; ++i) {
            result[i][j] = matrix[i - row_count][j];
        }
    }
    return result;
}

// Kiszámolja a mátrix transzponáltját.
template<typename T>
Matrix<T> Matrix<T>::transposition() const {
    Matrix<T> result(column_count, row_count);
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < column_count; ++j) {
            result[j][i] = block[i][j];
        }
    }
    return result;
}

// Elvégzi a mátrixon a Gauss-Jordan eliminációt.
template<typename T>
Matrix<T> Matrix<T>::GaussJordanElimination() const {
    Matrix<T> result(*this);
    int finished = 0;
    int zero_row = row_count;

    for (int j = 0; j < column_count && finished < zero_row; ++j) {
        bool zero_column = true;
        for (int i = finished; i < row_count && finished < zero_row; ++i) {
            if (result[i].isVectorZero()) {
                result.swapRows(--zero_row, i);
            }
            if (result[i][j] != 0) {
                result.swapRows(finished, i);
                zero_column = false;
                i += row_count;
            }
        }

        if (finished >= zero_row) break;

        result[finished] = result[finished] / result[finished][j];
        for (int i = 0; i < row_count; ++i) {
            if (i != finished) {
                result[i] = result[i] - result[i][j] * result[finished];
            }
        }
        if (!zero_column) finished++;
    }

    return result;
}

// Kiírja a mátrix tartalmát egy fájlba.
template<typename T>
void Matrix<T>::writeToFile(const std::string& filename) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw MatrixException("Hiba a fájl megnyitásakor: " + filename);
    }

    // Mátrix dimenzióinak írása a fájlba
    file.write(reinterpret_cast<const char*>(&row_count), sizeof(unsigned));
    file.write(reinterpret_cast<const char*>(&column_count), sizeof(unsigned));

    // Mátrix elemeinek írása a fájlba
    for (unsigned i = 0; i < row_count; ++i) {
        const Vektor<T>& row = block[i];
        file.write(reinterpret_cast<const char*>(row.getData()),
                   column_count * sizeof(T));
    }

    file.close();
}

// Beolvassa a mátrix tartalmát egy fájlból.
template<typename T>
void Matrix<T>::readFromFile(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw MatrixException("Hiba a fájl megnyitásakor: " + filename);
    }

    // Mátrix dimenzióinak beolvasása a fájlból
    file.read(reinterpret_cast<char*>(&row_count), sizeof(unsigned));
    file.read(reinterpret_cast<char*>(&column_count), sizeof(unsigned));

    // Mátrix foglalása és elemeinek beolvasása a fájlból
    block.resize(row_count);
    for (unsigned i = 0; i < row_count; ++i) {
        Vektor<T>& row = block[i];
        row.resize(column_count);
        file.read(reinterpret_cast<char*>(row.getData()), column_count * sizeof(T));
    }

    file.close();
}

// Túlterhelt indexelő operátorok, amelyek lehetővé teszik a mátrix elemeinek elérését.
template<typename T>
Vektor<T> Matrix<T>::operator[](int i) const {
    if (i < 0 || i >= row_count) {
        throw MatrixException("Érvénytelen sorindex.");
    }
    return block[i];
}

template<typename T>
Vektor<T>& Matrix<T>::operator[](int i) {
    if (i < 0 || i >= row_count) {
        throw MatrixException("Érvénytelen sorindex.");
    }
    return block[i];
}

// Mátrixműveletek támogatása operátorok túlterhelésével.
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& matrix) const {
    if (row_count != matrix.row_count || column_count != matrix.column_count) {
        throw MatrixException("A mátrixok dimenziói nem egyeznek a kivonás műveleténél.");
    }
    Matrix<T> result(row_count, column_count);
    result.block = block - matrix.block;
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrix) const {
    if (row_count != matrix.row_count || column_count != matrix.column_count) {
        throw MatrixException("A mátrixok dimenziói nem egyeznek az összeadás műveleténél.");
    }
    Matrix<T> result(row_count, column_count);
    result.block = block + matrix.block;
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrix) const {
    if (column_count != matrix.row_count) {
        throw MatrixException("A bal oldali mátrix oszlopainak száma nem egyezik a jobb oldali mátrix sorainak számával.");
    }
    Matrix<T> result(row_count, matrix.column_count);
    for (int i = 0; i < row_count; ++i) {
        for (int j = 0; j < matrix.column_count; ++j) {
            for (int k = 0; k < column_count; ++k) {
                result[i][j] += block[i][k] * matrix[k][j];
            }
        }
    }
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const {
    Matrix<T> result(row_count, column_count);
    for (int i = 0; i < row_count; ++i) {
        result[i] = block[i] * scalar;
    }
    return result;
}

// Barát függvény, ami lehetővé teszi a mátrix kiírását egy kimeneti adatfolyamba.
template<typename U>
std::ostream& operator<<(std::ostream& ostream, const Matrix<U>& matrix) {
    for (int i = 0; i < matrix.row_count; ++i) {
        ostream << "|";
        for (int j = 0; j < matrix.column_count; ++j) {
            ostream << matrix.block[i][j];
            if (j != matrix.column_count - 1) {
                ostream << " ";
            }
        }
        ostream << "|\n";
    }
    return ostream;
}


#endif // HAZI_MATRIX_MATRIX_H