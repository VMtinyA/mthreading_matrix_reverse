#include <QCoreApplication>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>

#include <QVector>

#define MAX_THREADS 2
#define LEN 128
#define ERROR_OPEN_FILE -3

using namespace std;

double det(const QVector<QVector<double> > &matrix, int size); // вычисляет определитель матрицы
void * AlgebraicComplement (void* thread_data);
bool MThReverse(QVector <QVector <double>> &matrix, int size);
void InitMatrixf(const char fname [LEN], QVector <QVector <double>> &matrix, const int n);
void ImportMatrixf(const char *fname, QVector <QVector <double>> &matrix, const int n);

//специальная структура для данных потока
typedef struct{
    int row_index; //номер обрабатываемой строки
    int col_index; // номер обрабатываемого элемента
    int size; //размер исходной матрицы
    //указатели на матрицы
    QVector <QVector <double>> copy;
    double adj; // возврат ответа
} pthrData;


int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    int n;
    printf("Enter n\n");
    scanf("%i", &n);

    QVector <QVector <double>> matrix;

    char fname [LEN];
    while (getc(stdin) != '\n');
    printf("Enter input filename: ");
    gets (fname);
    InitMatrixf(fname, matrix, n);

    MThReverse(matrix, n);

    printf("Enter output filename: ");
    gets (fname);
    ImportMatrixf(fname, matrix, n);

    exit(0);
    return a.exec();
}

void InitMatrixf(const char fname [LEN], QVector <QVector <double>> &matrix, const int n) // ввод из файла библиотекой С
{
    FILE *input;

    input = fopen (fname, "rb");
    if (input == NULL) exit(ERROR_OPEN_FILE);

    double temp = 0;
    matrix.resize(n);
    for (int i = 0; i < n; i++)
    {
        matrix [i].resize(n);
        for (int j = 0; j < n; j++)
        {
            fscanf(input, "%lf", &temp);
            matrix [i][j] = temp;
        }
    }

    if (fclose(input)) printf ("Closing Error");
}

void ImportMatrixf(const char *fname, QVector <QVector <double>> &matrix, const int n)
{
    FILE *output;

    output = fopen (fname, "w");
    if (output == NULL) exit(ERROR_OPEN_FILE);

    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++) fprintf(output, "%.3f ", matrix [i][j]);
        fprintf(output,"\n");
    }

    if (fclose(output)) printf ("Closing Error");
}

double det(const QVector <QVector <double>> &matrix, int size)
{
    //Функция, вычисляющая определитель матрицы
        // Принимает:
        // М - указатель на двумерный массив
        //  size - размерность матрицы
        // Возвращает:
        // double - число
        // Вычисление выполняется с помощью приведения к ступенчатому виду

    QVector <QVector <double>> M;
    M.resize(size);
    for (int i = 0; i < size; ++i)
    {
        M[i].resize(size);
        for (int j = 0; j < size; ++j)
            M[i][j] = matrix[i][j];
    }

    double s = 0;
    for (int k = 0; k < size; ++k)
    {
            // Если элемент на главной диагонали в исходной
            // строке - нуль, то ищем строку, где элемент
            // того же столбца не нулевой, и меняем строки
            // местами
         if (M[k][k] == 0)
         {
             bool changed = false;
             // Идём по строкам, расположенным ниже исходной
             for (int i = k + 1; i < size; ++i)
             {
                    // Если нашли строку, где в том же столбце
                    // имеется ненулевой элемент
                 if (M[i][k] != 0)
                 {
                        // Меняем найденную и исходную строки местами
                     swap(M[k], M[i]);
                     changed = true;
                     s++;
                     break;
                 }
             }
                // Если обмен строк произведён не был - матрица не может быть
                // обращена
             if (!changed)
                 return 0;
         }
         // Идём по строкам, которые расположены ниже исходной
         for (int i = k + 1; i < size; ++i)
         {
             // Запоминаем множитель
            double multi = M[i][k] / M[k][k];
             // Отнимаем от очередной строки исходную, умноженную
             // на сохранённый ранее множитель
            for (int j = 0; j < size; ++j)
                M[i][j]   -= multi * M[k][j];
         }
     }

     double ref = M[0][0];
     for (int i = 1; i < size; ++i)
        ref *= M[i][i];
     //delete *M;

     return ref * pow(-1, s);
}

void* AlgebraicComplement(void* thread_data)
{
    // Функция для нахождения алгебраического дополнения элемента (для многопоточности)
    // принимает и возвращает указатель на void
    // матрицу принимает из copy

    //получаем структуру с данными
    pthrData *data = (pthrData*) thread_data;
    // "выгрузка" аргументов
    int Msize = data->size - 1, row_index = data->row_index, col_index = data->col_index;

    // Создание и заполнение минора
    //double **Minor = new double *[Msize]();
    QVector <QVector <double>> Minor;
    int temp_row = 0, temp_col;
    Minor.resize(Msize);
    for (int Mi = 0; Mi < Msize; Mi++)
    {
        Minor[Mi].resize(Msize);
        if (temp_row == row_index) temp_row++;
        temp_col = 0;
        for (int Mj = 0; Mj < Msize; Mj++)
        {
            if (temp_col == col_index) temp_col++;
            Minor[Mi][Mj] = data->copy[temp_row][temp_col];
            temp_col++;
        }
        temp_row++;
    }

    // Запись ответа в adj
    data->adj = pow(-1, (row_index + col_index)) * det(Minor, Msize);
    //&thread_data = (void*) data;

    return NULL;
}

bool MThReverse(QVector <QVector <double>> &matrix, int size)
{
    //Функция, производящая обращение матрицы в несколько потоков
        // Принимает:
        //  matrix - указатель на массив массивов
        //  size - размерность матрицы
        // Возвращает:
        // true - успех; false - неудача.
        // Обращение выполняется с помощью присоединенной матрицы и det

    QVector <QVector <double>> copy;
        // Заполняем копию исходной матрицы
    copy.resize(size);
    for (int i = 0; i < size; ++i)
    {
        copy[i].resize(size);
        for (int j = 0; j < size; ++j)
            copy[i][j] = matrix[i][j];
    }

    double deter = det(copy, size);
    if (!deter) return false;

    // Транспонирование copy
    double temp;
    for (int i = 0; i < (size - 1); ++i)
        for (int j = i + 1; j < size; ++j)
        {
            temp = copy[j][i];
            copy[j][i] = copy[i][j];
            copy[i][j] = temp;
        }

    // Создание шаблона присоединенной матрицы
    // функция потока будет вести сюда запись полученного значения для каждого элемента
    QVector <QVector <double>> Adj;
    Adj.resize(size);
    for (int i = 0; i < size; ++i)
    {
        Adj[i].resize(size);
        for (int j = 0; j < size; ++j)
            Adj[i][j] = copy[i][j];
    }

    // Вычисление матрицы алгебраических дополнений в 2 потока (изменяет copy)

    pthread_t* threads = (pthread_t*) malloc(MAX_THREADS * sizeof(pthread_t));
    //pthrData *threadData = (pthrData*) malloc(MAX_THREADS * sizeof(pthrData));

    pthrData threadData[MAX_THREADS];
    int j, th_count;
    for (int i = 0; i < size; i++)
    {
        j = 0;
        while (j < size)
        {
            th_count = 0;
            while ( (th_count < MAX_THREADS) && (j != size) )
            {

                threadData[th_count].size = size;
                threadData[th_count].row_index = i;
                threadData[th_count].col_index = j;
                threadData[th_count].copy.resize(size);
                for (int x = 0; x < size; ++x)
                {
                    threadData[th_count].copy[x].resize(size);
                    for (int y = 0; y < size; ++y)
                        threadData[th_count].copy[x][y] = copy[x][y];
                }
                pthread_create (&(threads[th_count]), NULL, &AlgebraicComplement, &threadData[th_count]);
                th_count++;
                j++;
            }
            for (int k = th_count; k > 0; --k)
                pthread_join (threads[k - 1], NULL);
            for (int k = 0; k < th_count; ++k)
                Adj[threadData[k].row_index][threadData[k].col_index] = threadData[k].adj;
        }
    }
    //delete *copy;

    free(threads);

    // Итоговое вычисление обратной матрицы (делит каждый элемент copy на det исходной)
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
           matrix[i][j] = Adj[i][j] / deter;

    printf("Reversed\n");
    return true;
}
