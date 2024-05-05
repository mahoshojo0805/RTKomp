#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<iostream>
#include"omp.h"
#include"RTK_Struct.h"
#define T_ZERO 1e-10//判断为0的阈值
#define MATRIX_M 300//测试矩阵的行数和列数
#define THRESHOLD 0.001//迭代停止的精度阈值，单位为米
#define MAXCOUNT 10//迭代的最大次数
#pragma region 其余辅助函数
/**
 * @brief                   数组初始化为value
 * @param size
 * @param value
 * @param dst
 * @note
 */
void SetArray(const int size, const double value, double dst[])
{
    for (int i = 0; i < size; i++)
        dst[i] = value;
}

/**
 * @brief                   复制数组串行实现
 * @param[in]   size        数组大小
 * @param[in]   src         源数组
 * @param[out]  dst         目标数组
 * @note
 */
void CopyArray(const int size, const double src[], double dst[])
{
    for (int i = 0; i < size; i++)
        dst[i] = src[i];
}
/**
 * @brief                   计算两点之间的集合距离
 * @param X1
 * @param X2
 * @return                  距离
 * @note
 */
double GetDistance(double X1[], double X2[])
{
    return sqrt(pow(X1[0] - X2[0], 2) + pow(X1[1] - X2[1], 2) + pow(X1[2] - X2[2], 2));
}
#pragma endregion

#pragma region 文件操作
/**
 * @brief                   读取卫星位置
 * @param file
 * @param data
 * @return                  true：读取成功 false：读取失败或者已到文件末尾
 * @note
 */
bool ReadSatPos(FILE* file, Data* data)
{
    char c = 'E';
    data->satnum[0] = 0;
    data->satnum[1] = 0;
    bool first[2] = {true,true};//某系统的第一个卫星选为参考星
    while (true)
    {
        fscanf_s(file, "%c", &c);
        if (feof(file) != 0)return false;
        if (c == '#')
        {
            fscanf_s(file, "%d %lf", &data->time.week, &data->time.sow);
            fscanf_s(file, "%c", &c);
        }
        else if (c == 'E')
        {
            fscanf_s(file, "%c", &c);
            break;
        }
        else {
            int sys = (c == 'G' ? GPS : BDS);
            int idx = data->satnum[0] + data->satnum[1];
            data->satnum[sys]++;
            data->satpos[idx].sys = (c == 'G' ? GPS : BDS);
            if (first[data->satpos[idx].sys])
            {
                data->refsatidx[data->satpos[idx].sys] = idx;
                first[data->satpos[idx].sys] = false;
            }
            fscanf_s(file, "%d %lf %lf %lf %lf %lf %lf", &data->satpos[idx].prn,
                &data->satpos[idx].posrov[0], &data->satpos[idx].posrov[1], &data->satpos[idx].posrov[2],
                &data->satpos[idx].posbas[0], &data->satpos[idx].posbas[1], &data->satpos[idx].posbas[2]);
            fscanf_s(file, "%c", &c);
        }
    }
    return true;
}
/**
 * @brief                   读取单差观测值文件
 * @param file
 * @param data
 * @return                  true：文件读取成功 false:文件读取失败或者已到末尾
 * @note
 */
bool ReadObs(FILE* file, Data* data)
{
    char c='E';
    data->satnum[0] = 0;
    data->satnum[1] = 0;
    while (true)
    {
        fscanf_s(file, "%c", &c);
        if (feof(file) != 0)return false;
        if (c == '#')
        {
            fscanf_s(file, "%d %lf", &data->time.week, &data->time.sow);
            fscanf_s(file, "%c", &c);
        }
        else if (c == 'E')
        {
            fscanf_s(file, "%c", &c);
            break;
        }
        else {
            int sys = (c == 'G' ? GPS : BDS);
            int idx = data->satnum[0] + data->satnum[1];
            data->satnum[sys]++;
            data->obs[idx].sys = (c == 'G' ? GPS : BDS);
            fscanf_s(file, "%d %lf %lf %lf %lf", &data->obs[idx].prn, &data->obs[idx].P[0], &data->obs[idx].P[1],
                &data->obs[idx].L[0], &data->obs[idx].L[1]);
            fscanf_s(file, "%c", &c);
        }
    }
    return true;
}
/**
 * @brief                   将矩阵写到文件中
 * @param filename
 * @param m
 * @param n
 * @param mat
 * @note
 */
void WriteMatrix(const char* filename, const int m, const int n, const double mat[])
{
    FILE* file;//文件指针
    fopen_s(&file, filename, "w");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fprintf(file, "%10.3f\t", mat[i * n + j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

/**
 * @brief                   从文件中读取矩阵
 * @param filename
 * @param m
 * @param n
 * @param mat
 * @note
 */
void ReadMatrix(const char* filename, const int m, const int n, double mat[])
{
    FILE* file;
    fopen_s(&file, filename, "r");
    if (file == NULL)
    {
        printf("open failed\n");
        return;
    }
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            fscanf_s(file, "%lf", &mat[i * n + j]);
        }
    }
    fclose(file);
}

/**
 * @brief                   随机生成一个m*m的矩阵，并写入文件
 * @param filename
 * @param m
 * @note
 */
void RandMatrix(const char* filename, const int m, const int n)
{
    FILE* file;//文件指针
    fopen_s(&file, filename, "w");
    int min = -10, max = 10;
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double randnum = ((double)(rand() % (max * 1000 - min * 1000) + min * 1000) / 1000.0);
            fprintf(file, "%*.3f\t", 8, randnum);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}
#pragma endregion

#pragma region 矩阵相关函数

/**
 * @brief                   将矩阵打印到屏幕上
 * @param[in]   m           矩阵总行数
 * @param[in]   n           矩阵总列入
 * @param[in]   mat         矩阵
 * @note
 */
void Matrix_Print(const int m, const int n, const double mat[])
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%.4f ", mat[i * n + j]);
        }
        printf("\n");
    }
}

/**
 * @brief                   构造单位矩阵
 * @param[in]   m           矩阵的大小为m*m
 * @param[out]  i_mat       输出的单位矩阵
 * @note
 */
void Matrix_I(const int m, double i_mat[])
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < m; j++)
            i_mat[i * m + j] = i == j ? 1.0 : 0.0;
    }
}

/**
 * @brief                   矩阵两行交换
 * @param[in]       r       矩阵的总行数
 * @param[in]       c       矩阵的总列数
 * @param[in]       m       交换的两行之一
 * @param[in]       n       交换的两行之一
 * @param[in|out]   mat     矩阵
 * @note
 */
void Matrix_SwagRow(const int r, const int c, const int m, const int n, double mat[])
{
    if (m < 0 || n < 0 || m >= r || n >= r)return;
    for (int i = 0; i < c; i++)
    {
        double temp = mat[m * c + i];
        mat[m * c + i] = mat[n * c + i];
        mat[n * c + i] = temp;
    }
}

/**
 * @brief                   矩阵加法
 * @param m
 * @param n
 * @param A
 * @param B
 * @param output
 * @note
 */
void MatrixPlus(const int m, const int n, const double A[], const double B[], double output[])
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            output[i * n + j] = A[i * n + j] + B[i * n + j];
        }
    }
}

/**
 * @brief                   使用OPENMP矩阵求逆并行化
 * @param[in]   m           矩阵的大小为m*m
 * @param[in]   input_mat   输入的矩阵
 * @param[out]  inv_mat     输出，矩阵的逆
 * @return      true        求逆成功
 *              false       矩阵不可逆
 * @note
 */
bool Matrix_Inv_omp(const int m, const double input_mat[], double inv_mat[])
{
    int i, j, k;
    double* mat;
    double* I;//单位矩阵
    I = (double*)malloc(sizeof(double) * m * m);
    mat = (double*)malloc(sizeof(double) * m * m);//分配内存
    CopyArray(m * m, input_mat, mat);//复制输入的数组用于计算
    Matrix_I(m, I);
    //化为上三角矩阵
    for (i = 0; i < m; i++)
    {
        if (abs(mat[i * m + i]) < T_ZERO)//如果为0，与下面的非0行交换
        {
            for (j = i + 1; j < m; j++)
            {
                if (mat[j * m + i] == 0)
                    continue;
                Matrix_SwagRow(m, m, i, j, mat);
                Matrix_SwagRow(m, m, i, j, I);
                break;
            }
            if (abs(mat[i * m + i]) < T_ZERO)
                return false;//还是0，没找到非0行，不可逆
        }
        double temp = mat[i * m + i];
        for (j = 0; j < m; j++)//首元素化为1
        {
            mat[i * m + j] /= temp;
            I[i * m + j] /= temp;
        }
#pragma omp parallel for shared(i,m,mat,I) private(j,k)
        for (j = i + 1; j < m; j++)//对每一行减去scale*第i行，使首个非0元素变为0
        {
            double scale = mat[j * m + i];
            for (k = 0; k < m; k++)
            {
                mat[j * m + k] -= scale * mat[i * m + k];
                I[j * m + k] -= scale * I[i * m + k];
            }
        }
    }
    //化为下三角矩阵
    for (i = m - 1; i >= 0; i--)
    {
#pragma omp parallel for shared(i,m,mat,I) private(j,k)
        for (j = i - 1; j >= 0; j--)
        {
            double scale = mat[j * m + i];
            for (k = 0; k < m; k++)
            {
                mat[j * m + k] -= scale * mat[i * m + k];
                I[j * m + k] -= scale * I[i * m + k];
            }
        }
    }
    CopyArray(m * m, I, inv_mat);
    free(mat);//释放内存
    free(I);
    return true;
}
/**
 * @brief                   矩阵求逆串行实现
 * @param[in]   m           矩阵的大小为m*m
 * @param[in]   input_mat   输入的矩阵
 * @param[out]  inv_mat     输出，矩阵的逆
 * @return      true        求逆成功
 *              false       矩阵不可逆
 * @note
 */
bool Matrix_Inv(const int m, const double input_mat[], double inv_mat[])
{
    int i, j, k;
    double* mat;
    double* I;//单位矩阵
    I = (double*)malloc(sizeof(double) * m * m);
    mat = (double*)malloc(sizeof(double) * m * m);//分配内存
    CopyArray(m * m, input_mat, mat);//复制输入的数组用于计算
    Matrix_I(m, I);
    //化为上三角矩阵
    for (i = 0; i < m; i++)
    {
        if (abs(mat[i * m + i]) < T_ZERO)//如果为0，与下面的非0行交换
        {
            for (j = i + 1; j < m; j++)
            {
                if (mat[j * m + i] == 0)
                    continue;
                Matrix_SwagRow(m, m, i, j, mat);
                Matrix_SwagRow(m, m, i, j, I);
                break;
            }
            if (abs(mat[i * m + i]) < T_ZERO)
                return false;//还是0，没找到非0行，不可逆
        }
        double temp = mat[i * m + i];
        for (j = 0; j < m; j++)//首元素化为1
        {
            mat[i * m + j] /= temp;
            I[i * m + j] /= temp;
        }
        for (j = i + 1; j < m; j++)//对每一行减去scale*第i行，使首个非0元素变为0
        {
            double scale = mat[j * m + i];
            for (k = 0; k < m; k++)
            {
                mat[j * m + k] -= scale * mat[i * m + k];
                I[j * m + k] -= scale * I[i * m + k];
            }
        }
    }
    //化为下三角矩阵
    for (i = m - 1; i >= 0; i--)
    {
        for (j = i - 1; j >= 0; j--)
        {
            double scale = mat[j * m + i];
            for (k = 0; k < m; k++)
            {
                mat[j * m + k] -= scale * mat[i * m + k];
                I[j * m + k] -= scale * I[i * m + k];
            }
        }
    }
    CopyArray(m * m, I, inv_mat);
    free(mat);//释放内存
    free(I);
    return true;
}
/**
 * @brief                   矩阵转置
 * @param[in]   m           矩阵行数
 * @param[in]   n           矩阵列数
 * @param[in]   mat         输入的矩阵
 * @param[out]  output      输出的矩阵
 * @note
 */
void MatrixT(const int m, const int n, const double mat[],double output[])
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            output[j * m + i] = mat[i * n + j];
        }
    }
}
/**
 * @brief                   矩阵乘法并行实现
 * @param m
 * @param l
 * @param n
 * @param A
 * @param B
 * @param Tflag             是否转置的标记，2维字符数组，取值为"TT"/"TN"/"NT"/"NN"，T表示转置，N表示正常不转置
 * @param output
 * @note
 */
void MatrixMultiply_omp(const int m, const int l, const int n,
    const double A[], const double B[], const char Tflag[], double output[])
{
    int i, j, k;
    double* a = (double*)malloc(sizeof(double) * m * l);
    double* b = (double*)malloc(sizeof(double) * l * n);
    Tflag[0] == 'T' ? MatrixT(l, m, A, a) : CopyArray(m * l, A, a);
    Tflag[1] == 'T' ? MatrixT(n, l, B, b) : CopyArray(l * n, B, b);
#pragma omp parallel for shared(m,n,l,a,b,output) private(i,j,k)
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            output[i * n + j] = 0;
            for (k = 0; k < l; k++)
            {
                output[i * n + j] += a[i * l + k] * b[k * n + j];
            }
        }
    }
    free(a);
    free(b);
}
/**
 * @brief                   矩阵连乘A*B*C
 * @param m
 * @param l1
 * @param l2
 * @param n
 * @param A
 * @param B
 * @param C
 * @param Tflag
 * @param output
 * @note
 */
void MatrixMultiply_omp(const int m, const int l1, const int l2, const int n,
    const double A[], const double B[], const double C[], const char Tflag[], double output[])
{
    double* temp = (double*)malloc(sizeof(double) * m * l2);
    MatrixMultiply_omp(m, l1, l2, A, B, Tflag, temp);
    char temp_flag[2] = { 'N',Tflag[2] };
    MatrixMultiply_omp(m, l2, n, temp, C, temp_flag, output);
    free(temp);
}
/**
 * @brief                   矩阵乘法A*B，A/AT为m*l，B/BT为l*n
 * @param[in]   m                 
 * @param[in]   l
 * @param[in]   n
 * @param[in]   A           
 * @param[in]   B
 * @param[in]   Tflag       是否转置的标记，2维字符数组，取值为"TT"/"TN"/"NT"/"NN"，T表示转置，N表示正常不转置
 * @param[out]  output      相乘结果
 * @note
 */
void MatrixMultiply(const int m, const int l, const int n, 
    const double A[], const double B[], const char Tflag[], double output[])
{ 
    double* a = (double*)malloc(sizeof(double) * m * l);
    double* b = (double*)malloc(sizeof(double) * l * n);
    Tflag[0] == 'T' ? MatrixT(l, m, A, a) : CopyArray(m * l, A, a);
    Tflag[1] == 'T' ? MatrixT(n, l, B, b) : CopyArray(l * n, B, b);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            output[i * n + j] = 0;
            for (int k = 0; k < l; k++)
            {
                output[i * n + j] += a[i * l + k] * b[k * n + j];
            }
        }
    }
    free(a);
    free(b);
}
/**
 * @brief                   三个矩阵连乘的串行实现
 * @param m
 * @param l1
 * @param l2
 * @param n
 * @param A
 * @param B
 * @param C
 * @param Tflag
 * @param output
 * @note
 */
void MatrixMultiply(const int m, const int l1, const int l2, const int n,
    const double A[], const double B[], const double C[], const char Tflag[], double output[])
{
    double* temp = (double*)malloc(sizeof(double) * m * l2);
    MatrixMultiply(m, l1, l2, A, B, Tflag, temp);
    char temp_flag[2] = { 'N',Tflag[2] };
    MatrixMultiply(m, l2, n, temp, C, temp_flag, output);
    free(temp);
}
#pragma endregion

/**
 * @brief               读取文件并进行RTK
 * @param[in] omp       是否使用OpenMP加速
 * @param[in] filepath  保存解算结果的路径
 * @note
 */
void RTK(const bool omp, const char filepath[])
{
    Data data;
    FILE* obsfile;
    FILE* satfile;
    FILE* resfile;
    fopen_s(&obsfile, "Data\\SDOBSDATA.txt", "r");
    fopen_s(&satfile, "Data\\SATDATA.txt", "r");
    fopen_s(&resfile, filepath, "w");
    double* B = (double*)malloc(sizeof(double) * MAXM * MAXN);
    double* L = (double*)malloc(sizeof(double) * MAXM * 1);
    double* X = (double*)malloc(sizeof(double) * MAXN * 1);
    double* dX = (double*)malloc(sizeof(double) * MAXN * 1);
    double* tX = (double*)malloc(sizeof(double) * MAXN * 1);
    double* P = (double*)malloc(sizeof(double) * MAXM * MAXM);
    double* N = (double*)malloc(sizeof(double) * MAXN * MAXN);
    double* NInv = (double*)malloc(sizeof(double) * MAXN * MAXN);;
    double* W = (double*)malloc(sizeof(double) * MAXN * 1);
    double WL[2][2] = { {WL1_GPS,WL2_GPS},{WL1_BDS,WL3_BDS} };//波长
    double delta[2] = { 0.3,0.01 };//伪距和相位的误差
    double BasPos[3] = { -2267804.526310, 5009342.372351, 3220991.863215 };//基站精确坐标
    while (true)
    {
        //读取文件
        if (!ReadObs(obsfile, &data) || !ReadSatPos(satfile, &data))
            break;
        int ddnum = (data.satnum[0] + data.satnum[1] - 2);
        int m = ddnum * 4;
        int n = ddnum * 2 + 3;
        int Idx;
        int count;
        //初始化
        Idx = 0;
        SetArray(MAXN * 1, 0, X); //memcpy(X, BasPos, sizeof(double) * 3);
        SetArray(MAXN * 1, 0, dX);
        count = 0;
        //最小二乘
        do {
            SetArray(MAXM * MAXN, 0, B);
            SetArray(MAXM * MAXM, 0, P);
            SetArray(MAXM * 1, 0, L);
            Idx = 0;
            for (int i = 0; i < (data.satnum[0] + data.satnum[1]); i++)
            {
                int sysIdx = data.obs[i].sys;
                if (i == data.refsatidx[sysIdx])
                    continue;//参考星
                int refidx = data.refsatidx[sysIdx];
                double roR2 = GetDistance(X, data.satpos[i].posrov);//流动站到当前卫星的几何距离
                double roR1 = GetDistance(X, data.satpos[refidx].posrov);//流动站到参考星的几何距离
                double roB2 = GetDistance(BasPos, data.satpos[i].posbas);//基站到当前卫星的几何距离
                double roB1 = GetDistance(BasPos, data.satpos[refidx].posbas);//基站到参考星的几何距离
                double DDro = roR2 - roB2 - roR1 + roB1;//双差

                double lDD = (X[0] - data.satpos[i].posrov[0]) / roR2
                    - (X[0] - data.satpos[refidx].posrov[0]) / roR1;
                double mDD = (X[1] - data.satpos[i].posrov[1]) / roR2
                    - (X[1] - data.satpos[refidx].posrov[1]) / roR1;
                double nDD = (X[2] - data.satpos[i].posrov[2]) / roR2
                    - (X[2] - data.satpos[refidx].posrov[2]) / roR1;

                for (int j = 0; j < 2; j++)
                {
                    /*伪距*/
                    B[Idx * n + 0] = lDD;
                    B[Idx * n + 1] = mDD;
                    B[Idx * n + 2] = nDD;

                    L[Idx] = data.obs[i].P[j] - data.obs[refidx].P[j] - DDro;
                    if (sysIdx == 0)
                        for (int k = Idx % 4; k < (data.satnum[0] - 1) * 4; k += 4)
                            P[Idx * m + k] = k == Idx ? (double)(data.satnum[0] - 1) / 2 / delta[0] / (data.satnum[0])
                            : -1.0 / 2 / delta[0] / (data.satnum[0]);
                    else if (sysIdx == 1)
                        for (int k = Idx % 4 + (data.satnum[0] - 1) * 4; k < m; k += 4)
                            P[Idx * m + k] = k == Idx ? (double)(data.satnum[1] - 1) / 2 / delta[0] / (data.satnum[1])
                            : -1.0 / 2 / delta[0] / (data.satnum[1]);
                    Idx++;
                    /*相位*/
                    B[Idx * n + 0] = lDD;
                    B[Idx * n + 1] = mDD;
                    B[Idx * n + 2] = nDD;
                    double lamda = WL[sysIdx][j];//sysIdx系统的j频率的波长
                    B[Idx * n + 3 + Idx / 2] = lamda;
                    L[Idx] = (data.obs[i].L[j] - data.obs[refidx].L[j]) * lamda - DDro - lamda * X[3 + Idx / 2];
                    if (sysIdx == 0)
                        for (int k = Idx % 4; k < (data.satnum[0] - 1) * 4; k += 4)
                            P[Idx * m + k] = k == Idx ? (double)(data.satnum[0] - 1) / 2 / delta[1] / (data.satnum[0])
                            : -1.0 / 2 / delta[1] / (data.satnum[0]);
                    else if (sysIdx == 1)
                        for (int k = Idx % 4 + (data.satnum[0] - 1) * 4; k < m; k += 4)
                            P[Idx * m + k] = k == Idx ? (double)(data.satnum[1] - 1) / 2 / delta[1] / (data.satnum[1])
                            : -1.0 / 2 / delta[1] / (data.satnum[1]);
                    Idx++;
                }
            }
            omp ? MatrixMultiply_omp(n, m, m, n, B, P, B, "TNN", N) : MatrixMultiply(n, m, m, n, B, P, B, "TNN", N);
            bool inv_result = omp ? Matrix_Inv_omp(n, N, NInv) : Matrix_Inv(n, N, NInv);
            if (!inv_result)
            {
                printf("求逆失败\n");
                break;
            }
            omp ? MatrixMultiply_omp(n, m, m, 1, B, P, L, "TNN", W) : MatrixMultiply(n, m, m, 1, B, P, L, "TNN", W);
            omp ? MatrixMultiply_omp(n, n, 1, NInv, W, "NN", dX) : MatrixMultiply(n, n, 1, NInv, W, "NN", dX);
            MatrixPlus(n, 1, X, dX, tX);
            CopyArray(n * 1, tX, X);
            count++;
        } while (sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]) >= THRESHOLD && count <= MAXCOUNT);
        if (sqrt(dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]) < THRESHOLD && count <= 10)
            fprintf_s(resfile, "%8d%16.3lf%16.3lf%16.3lf%16.3lf\n", data.time.week, data.time.sow, X[0], X[1], X[2]);
    }
    fclose(obsfile);
    fclose(satfile);
    fclose(resfile);
    free(B); free(L); free(X); free(dX); free(tX);
    free(P); free(N); free(NInv); free(W);
}

int main()
{
    srand(time(NULL));
    omp_set_num_threads(8);//设置使用的线程数量
    time_t t1, t2;
    double dt;//用于计时，计算加速比

#pragma region TestForInv
    //double* mat;
    //double* inv;
    //mat = (double*)malloc(sizeof(double) * MATRIX_M * MATRIX_M);
    //inv = (double*)malloc(sizeof(double) * MATRIX_M * MATRIX_M);
    //RandMatrix("Data\\InputMatrix.txt", MATRIX_M, MATRIX_M);//随机生成矩阵

    //ReadMatrix("Data\\InputMatrix.txt", MATRIX_M, MATRIX_M, mat);//读取文件
    //for (int i = 0; i < 300; i++) {
    //    t1 = clock();//使用OPENMP
    //    Matrix_Inv_omp(MATRIX_M, mat, inv);
    //    t2 = clock();
    //    dt = t2 - t1;
    //    printf("With Openmp costs %lf ms\n", dt);
    //    WriteMatrix("Data\\OMP_Result.txt", MATRIX_M, MATRIX_M, inv);

    //    t1 = clock();//不使用OPENMP
    //    Matrix_Inv(MATRIX_M, mat, inv);
    //    t2 = clock();
    //    dt = t2 - t1;
    //    printf("Without Openmp costs %lf ms\n", dt);
    //    WriteMatrix("Data\\Result.txt", MATRIX_M, MATRIX_M, inv);
    //}
#pragma endregion

#pragma region TestForMultiply
    //double* a = (double*)malloc(sizeof(double) * MAXM * MAXN);
    //double* b = (double*)malloc(sizeof(double) * MAXM * MAXM);
    //double* c = (double*)malloc(sizeof(double) * MAXM * MAXM);
    //RandMatrix("Data\\Mul_A.txt", MAXM, MAXN);//随机生成矩阵
    //RandMatrix("Data\\Mul_B.txt", MAXM, MAXM);
    ////printf("%d*%d*%d=%d", sizeof(double), MAXM, MAXN,sizeof(double)*MAXN);
    //ReadMatrix("Data\\Mul_A.txt", MAXM, MAXN, a);//读取文件
    //ReadMatrix("Data\\Mul_B.txt", MAXM, MAXM, b);

    //t1 = clock();
    ////MatrixMultiply_omp(MATRIX_M, MATRIX_M, MATRIX_M, a, b, "NN", c);
    //MatrixMultiply_omp(MAXN, MAXM, MAXM, MAXN, a, b, a, "TNN", c);
    //t2 = clock();
    //dt = t2 - t1;
    //printf("With Openmp costs %lf ms\n", dt);

    ////WriteMatrix("Data\\Mul_Result_omp.txt", MATRIX_M, MATRIX_M, c);

    //t1 = clock();
    ////MatrixMultiply(MATRIX_M, MATRIX_M, MATRIX_M, a, b, "NN", c);
    //MatrixMultiply(MAXN, MAXM, MAXM, MAXN, a, b, a, "TNN", c);
    //t2 = clock();
    //dt = t2 - t1;
    //printf("Without Openmp costs %lf ms\n", dt);

    //double testa[6] = { 3,2,5,1,3,6 };
    //double testb[9] = { 5,1,2,3,6,4,8,7,6 };
    //double testc[4];
    //MatrixMultiply_omp(2, 3, 3, 2, testa, testb, testa, "TNN", testc);
    //Matrix_Print(3, 2, testa);
    //Matrix_Print(3, 3, testb);
    //Matrix_Print(2, 2, testc);
#pragma endregion

#pragma region LSQ
    t1 = clock();
    RTK(true,"Data\\POS_OMP.txt");
    t2 = clock();
    dt = t2 - t1;
    printf("With Openmp cost %lf ms\n", dt);

    t1 = clock();
    RTK(false,"Data\\POS.txt");
    t2 = clock();
    dt = t2 - t1;
    printf("Without Openmp cost %lf ms\n", dt);
#pragma endregion

    return 0;
}
