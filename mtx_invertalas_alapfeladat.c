#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
Nagy Péter, m07ilf
2017.03.10.
*/

const double threshold=0.000000000001;

struct StrMtx{
    int size;
    double* Mtx;
};

struct Element{
    int RowPlace;
    int ColPlace;
    double value;
};

struct StrElimMtx{
    struct StrMtx Left;
    struct StrMtx Right;
};

void Usage(){



    printf("Hasznalat (invertalas): ./a.out matrixMeret matrixSrc \n");
    printf("vagy (egyenletrendszer): ./a.out matrixMeret matrixSrc vektorSrc\n");
    exit(0);
}
double *ReadMatrix(char* filename, int size){
    int i,j;
    FILE *inMatrixStream;
    inMatrixStream=fopen(filename, "r");
    double *matrix;
    matrix=(double *)malloc(size*size*sizeof(double));
    for(i=0; i<size;i++){
        for(j=0;j<size;j++){
            int dataRead=fscanf(inMatrixStream, "%lf", &matrix[i*size+j]);
            if(dataRead==0){
                printf("Hiba a beolvasasnal. \n");
            }
        }
    }
    return matrix;
}

void PrintMatrix(struct StrMtx inputmatrix){
    int i,j;
    for(i=0;i<inputmatrix.size;i++){
        for(j=0;j<inputmatrix.size;j++){
            printf("%f ", inputmatrix.Mtx[inputmatrix.size*i+j]);
        }
        printf("\n");
    }
}

double* Szorzas(double *matrix, double* vector, int rows, int cols){
    int i, j;
    double *productvec=malloc(rows*sizeof(double));
    for(i=0;i<rows;i++){
        double sum=0.0;
        for(j=0;j<cols;j++){
            sum=sum+matrix[i*cols+j]*vector[j];
        }
        productvec[i]=sum;
    }
    return productvec;
}

struct StrMtx IdN(int n){
    struct StrMtx idn;
    double *idnmatrix;
    int i;
    int j;
    idnmatrix= (double*) malloc(sizeof(double)*n*n);
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i==j){
                idnmatrix[i*n+j]=1;
            }
            else{
                idnmatrix[i*n+j]=0;
            }
        }
    }
    idn.size=n;
    idn.Mtx=idnmatrix;
    return idn;
}

void SubtractMultipleRows(struct StrMtx inputmatrix,double Coefficient, int lineNumber1, int lineNumber2){
    int i=0;
    for(i=0;i<inputmatrix.size; i++){
        inputmatrix.Mtx[lineNumber1*inputmatrix.size+i]-=Coefficient*inputmatrix.Mtx[lineNumber2*inputmatrix.size+i];
    }
}

void SwapRows(struct StrMtx inputmatrix, int lineNumber1, int lineNumber2){
    int i=0;
    double temp1=0.0;
    double temp2=0.0;
    for(i=0;i<inputmatrix.size; i++){
        temp1=inputmatrix.Mtx[lineNumber1*inputmatrix.size+i];
        temp2=inputmatrix.Mtx[lineNumber2*inputmatrix.size+i];
        inputmatrix.Mtx[lineNumber1*inputmatrix.size+i]=temp2;
        inputmatrix.Mtx[lineNumber2*inputmatrix.size+i]=temp1;
    }

}

void SwapCols(struct StrMtx inputmatrix, int colNumber1, int colNumber2){
    int i=0;
    if(colNumber1==colNumber2){
        return;
    }
    double temp1=0.0;
    double temp2=0.0;
    for(i=0;i<inputmatrix.size; i++){
        temp1=inputmatrix.Mtx[i*inputmatrix.size+colNumber1];
        temp2=inputmatrix.Mtx[i*inputmatrix.size+colNumber2];
        inputmatrix.Mtx[i*inputmatrix.size+colNumber2]=temp1;
        inputmatrix.Mtx[i*inputmatrix.size+colNumber1]=temp2;
    }

}

struct Element GetMaxInRow(struct StrMtx inputmatrix, int lineNumber, int startpoint){
    int i;
    struct Element result;
    result.value=inputmatrix.Mtx[lineNumber*inputmatrix.size+startpoint];
    result.RowPlace=lineNumber;
    result.ColPlace=startpoint;

    for(i=startpoint; i<inputmatrix.size; i++){
        if(fabs(inputmatrix.Mtx[inputmatrix.size*lineNumber+i])> fabs(result.value)){
            result.value=inputmatrix.Mtx[inputmatrix.size*lineNumber+i];
            result.ColPlace=i;
        }
    }
    return result;
}

struct Element GetMaxInSubMtx(struct StrMtx inputmatrix, int lineNumber, int colNumber){
    int i;
    struct Element result;
    struct Element temp;
    result.value=inputmatrix.Mtx[lineNumber*inputmatrix.size+colNumber];
    result.RowPlace=lineNumber;
    result.ColPlace=colNumber;
    temp=result;
    for(i=lineNumber; i<inputmatrix.size; i++){
        temp=GetMaxInRow(inputmatrix, i, colNumber);
        if(fabs(temp.value)>fabs(result.value)){
            result=temp;
        }
    }
    return result;
}

void MultiplyRow(struct StrMtx inputmatrix, int lineNumber, double Coefficient){
    int i;
    for(i=0;i<inputmatrix.size;i++){
        inputmatrix.Mtx[lineNumber*inputmatrix.size+i]=inputmatrix.Mtx[lineNumber*inputmatrix.size+i]*Coefficient;
    }
}

void Norm(struct StrElimMtx matrix){
    int j;
    for(j=0;j<matrix.Left.size;j++){
        double MaxElementValue;
        MaxElementValue=GetMaxInRow(matrix.Left, j, 0).value;
        if(MaxElementValue<threshold){
            printf("A matrix szingularis/rosszul kondicionalt \n");
            exit(0);
        }
        MultiplyRow(matrix.Left, j, 1./(MaxElementValue));
        MultiplyRow(matrix.Right, j, 1./(MaxElementValue));
    }
}

void Pivot(struct StrElimMtx matrix, int rownumber, int* ColPermutation){
    int n=matrix.Left.size;
    if(matrix.Left.Mtx[rownumber*n+rownumber]<threshold){
        struct Element MaxValue=GetMaxInSubMtx(matrix.Left, rownumber, rownumber);
        if(fabs(MaxValue.value)>=threshold){
            SwapRows(matrix.Left, rownumber, MaxValue.RowPlace);
            SwapRows(matrix.Right, rownumber, MaxValue.RowPlace);
            SwapCols(matrix.Left, rownumber, MaxValue.ColPlace);
            SwapCols(matrix.Right, rownumber, MaxValue.ColPlace);
            ColPermutation[rownumber]=MaxValue.ColPlace;
            }
        else{
            printf("A matrix szingularis \n");
            exit(0);
        }
    }
}

void Invert(struct StrElimMtx matrix){
    int n=matrix.Left.size;
    int k;
    int* permutation;
    permutation=(int *) malloc(n*sizeof(int));
    int d;
    for(d=0;d<n;d++){
        permutation[d]=d;
    }
    Norm(matrix);
    double diag;
    for(k=0;k<n;k++){
        if(fabs(matrix.Left.Mtx[k*n+k])<threshold){
            Pivot(matrix, k, permutation);
        }
        diag=matrix.Left.Mtx[k*n+k];
        int i;

        struct Element Maximum=GetMaxInRow(matrix.Left, k,0);
        if(fabs(Maximum.value)<threshold){
            printf("A matrix szingularis \n");
            exit(0);
        }
        for(i=k+1;i<n;i++){
            double ratio=matrix.Left.Mtx[i*n+k]/diag;
            SubtractMultipleRows(matrix.Left,ratio, i, k);
            SubtractMultipleRows(matrix.Right,ratio, i, k);
        }
        MultiplyRow(matrix.Left, k, 1./diag);
        MultiplyRow(matrix.Right, k, 1./diag);
    }

    for(k=n-1;k>=0;k--){
        int i;
        for(i=k-1;i>=0;i--){
            double ratio=matrix.Left.Mtx[i*n+k];
            SubtractMultipleRows(matrix.Left, ratio, i, k);
            SubtractMultipleRows(matrix.Right, ratio, i, k);
        }
    }
    for(d=0;d<n;d++){
        if(permutation[d]!=d){
            SwapCols(matrix.Left, permutation[d], d);
            SwapCols(matrix.Right, permutation[d], d);
            SwapRows(matrix.Left, permutation[d], d);
            SwapRows(matrix.Right, permutation[d], d);
        }
    }
}

int main(int argc , char *argv[]){
    if(argc!=3 && argc!=4){
        printf("%d", argc);

        Usage();
    }
    struct StrMtx argmatrix;
    argmatrix.size=atoi(argv[1]);
    argmatrix.Mtx=ReadMatrix(argv[2], atoi(argv[1]));

    struct StrMtx idn=IdN(argmatrix.size);
    struct StrElimMtx result;
    result.Left.size=argmatrix.size;
    result.Right.size=argmatrix.size;
    result.Left=argmatrix;
    result.Right=idn;
    printf("Beolvasott Matrix: \n");
    PrintMatrix(result.Left);
    Invert(result);
    if(argc==3){
        printf("\n Inverz: \n");
        PrintMatrix(result.Right);
    }

    if(argc==4){
        int i;
        int n=result.Left.size;
        double *InVector;
        FILE *InVectorStream;
        InVectorStream=fopen(argv[3], "r");
        InVector=(double *)malloc(n*sizeof(double));
        for(i=0;i<n;i++){
            fscanf(InVectorStream, "%lf", &InVector[i]);
        }
        double *ResVector;
        ResVector=Szorzas(result.Right.Mtx, InVector, n, n);
        printf("\n Megoldas: \n");
        for(i=0;i<n;i++){
            printf("%lf\n", ResVector[i]);
        }
    }
    return 0;
}
