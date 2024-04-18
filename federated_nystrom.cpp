#include<iostream>
#include<fstream>
#include<cctype>
#include<iomanip>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

MatrixXd computePositiveDifference(const MatrixXd& matrix1, const MatrixXd& matrix2) {

    MatrixXd difference = matrix1 - matrix2;


    difference = difference.array().abs();

    return difference;
}
MatrixXd calculatePercentageErrorMatrix(const MatrixXd& matrix1, const MatrixXd& matrix2) {
    int rows = matrix1.rows();
    int cols = matrix1.cols();

    MatrixXd percentageErrorMatrix(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double error = std::abs((matrix1(i, j) - matrix2(i, j)) / matrix1(i, j)) * 100.0;
            percentageErrorMatrix(i, j) = error;
        }
    }

    return percentageErrorMatrix;
}

double calculateAveragePercentageError(const MatrixXd& matrix1, const MatrixXd& matrix2) {
    double totalError = 0.0;
    int numElements = matrix1.size();

    for (int i = 0; i < matrix1.rows(); ++i) {
        for (int j = 0; j < matrix1.cols(); ++j) {
            double error = abs((matrix1(i, j) - matrix2(i, j)) / matrix1(i, j)) * 100.0;
            totalError += error;
        }
    }

    return totalError / numElements;
}
double calculateFrobeniusNormError(const MatrixXd& matrix1, const MatrixXd& matrix2) {
    MatrixXd differenceMatrix = matrix1 - matrix2; // Calculate the difference matrix
    double frobeniusNorm = differenceMatrix.norm(); // Calculate the Frobenius norm
    return frobeniusNorm;
}
int main() {
    label:
    int ch;
    //int num;
    //do
    //{
        system("CLS");
        cout << "\n\n\t\t\t\t======================\n";
        cout << "\t\t\t\tFEDERATED NYSTROM COMPLETION";
        cout << "\n\t\t\t\t======================\n";
        cout << "\t\t\t\t ::MAIN MENU::\n";
        cout << "\n\t\t\t\t1. NYSTROM APPROXIMATION";
        cout << "\n\t\t\t\t2. FEDERATED NYSTROM APPROACH USING FIXED CONSOLIDATED DATA";
        cout << "\n\t\t\t\t3. FEDERATED NYSTROM USING DIFFERENT KERNELS";
        cout << "\n\t\t\t\t4. EXIT";

        cout << "\n\n\t\t\t\tSelect Your Option (1-4): ";
        cin >> ch;


        if (ch == 1) {
            system("CLS");
            // Define the original matrix dimensions and rank
            int rank = 4;
            int size = 32;
            int rs = size - rank;

            // File handling to read matrix values from a file
            ifstream inputFile5("matrix_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile5.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

           // inputFile5 >> rank >> size;
            rs = size - rank;

            MatrixXd kernel(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile5 >> kernel(i, j);
                }
            }

            inputFile5.close();


            MatrixXd A = kernel.block(0, 0, rank, rank); // submatrix a (2x2)
            MatrixXd B = kernel.block(0, rank, rank, rs); // submatrix b (2x14)
            MatrixXd Qorg = kernel.block(rank, rank, rs, rs); // submatrix d (14x14)

            // Display submatrices a, b, c, and d
            cout << "Submatrix a :\n" << A << "\n\n";
            cout << "Submatrix b :\n" << B << "\n\n";
            cout << "Submatrix q :\n" << Qorg << "\n\n";
            MatrixXd Ain(rank, rank);
            MatrixXd Bt(rs, rank);
            MatrixXd Q(rs, rs);

            Ain = A.inverse();
            Bt = B.transpose();
            Q = Bt * Ain * B;

            cout << "final" << endl << Q << endl;
            MatrixXd positiveDifference = computePositiveDifference(Qorg, Q);
            cout << "Difference Matrix with Positive Values:\n" << positiveDifference << "\n";
            double avgPercentageError = calculateAveragePercentageError(Qorg, Q);
            cout << "Average Percentage Error: " << avgPercentageError << "%\n";
            double frobeniusNormError = calculateFrobeniusNormError(Qorg, Q);
            cout << "Frobenius Norm Error: " << frobeniusNormError << "\n";
            system("pause");

        }
        if (ch == 2) {
            system("CLS");
            // Define the original matrix dimensions and rank
            int rank = 4;
            int size = 32;
            int rs = size - rank;

            // File handling to read matrix values from a file
            ifstream inputFile0("matrix_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile0.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

            //inputFile0 >> rank >> size;
            rs = size - rank;

            MatrixXd kernel(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile0 >> kernel(i, j);
                }
            }

            inputFile0.close();

            // Display the kernel matrix
           // cout << "Kernel Matrix:\n" << kernel << "\n\n";

            // Define submatrices a, b, c, and d
            MatrixXd A = kernel.block(0, 0, rank, rank); // submatrix a (2x2)
            MatrixXd B = kernel.block(0, rank, rank, rs); // submatrix b (2x14)
            MatrixXd Qorg = kernel.block(rank, rank, rs, rs); // submatrix d (14x14)

            // Display submatrices a, b, c, and d
            /*cout << "Submatrix a :\n" << A << "\n\n";
            cout << "Submatrix b :\n" << B << "\n\n";
            cout << "Submatrix q :\n" << Qorg << "\n\n";*/
            MatrixXd Ain(rank, rank);
            MatrixXd Bt(rs, rank);
            MatrixXd Q(rs, rs);

            const int n = 3;
            // cout << "enter the no. of clients" << endl;
             //cin >> n;
            int row[n];
            int col[n];
            for (int i = 0; i < n; i++)
            {
                cout << "enter dimensions of client " << i + 1 << endl;
                cin >> row[i];
                if (row[i] > rank)
                {
                    cout << "invalid dimensions";
                    return 0;
                }
                cin >> col[i];
                if (col[i] > rs)
                {
                    cout << "invalid dimensions";
                    return 0;
                }
            }
            array<MatrixXd, n> subB;
            array<MatrixXd, n> subA;
            array<MatrixXd, n> subAin;
            array<MatrixXd, n> subBt;
            array<MatrixXd, n> subQ;
            array<MatrixXd, n> subQor;

            for (int i = 0; i < n; i++)
            {
                subB[i] = B.block(0, 0, row[i], col[i]);
                subA[i] = A.block(0, 0, row[i], row[i]);
                subAin[i] = subA[i].inverse();
                subBt[i] = subB[i].transpose();
                subQ[i] = subBt[i] * subAin[i] * subB[i];
                subQor[i] = Qorg.block(0, 0, col[i], col[i]);
                cout << "final of " << i + 1 << endl << subQ[i] << endl;
                MatrixXd positiveDifference = computePositiveDifference(subQor[i], subQ[i]);
                cout << "Difference Matrix of " << i + 1 << " client with Positive Values : \n" << positiveDifference << "\n";
                double avgPercentageError = calculateAveragePercentageError(subQor[i], subQ[i]);
                MatrixXd percentageErrorMatrix = calculatePercentageErrorMatrix(subQor[i], subQ[i]);
                cout << "Percentage Error Matrix of " << i + 1 << " client with :\n" << percentageErrorMatrix << "\n\n";
                cout << "Average Percentage Error of " << i + 1 << " client with :\n " << avgPercentageError << "%\n";
                double frobeniusNormError = calculateFrobeniusNormError(subQor[i], subQ[i]);
                cout << "Frobenius Norm Error: " << frobeniusNormError << "\n";
            }

            int maxCol = 0;
            for (int i = 1; i < n; ++i) {
                if (col[i] > col[maxCol]) {
                    maxCol = i;
                }
            }
            //cout << maxCol;
            MatrixXd Qag(col[maxCol], col[maxCol]);
            Qag = subQ[maxCol];

            int** nagg = new int* [col[maxCol]];
            for (int i = 0; i < col[maxCol]; ++i) {
                nagg[i] = new int[col[maxCol]];
            }

            for (int i = 0; i < col[maxCol]; ++i) {
                for (int j = 0; j < col[maxCol]; ++j) {
                    nagg[i][j] = 1;
                }
            }

            /*for (int i = 0; i < col[maxCol]; ++i) {
                for (int j = 0; j < col[maxCol]; ++j) {
                    cout<<nagg[i][j] <<", ";
                }
                cout << endl;
            }*/

            for (int x = 0; x < n; x++)
            {
                for (int i = 0; i < col[x]; i++)
                {
                    for (int j = 0; j < col[x]; j++)
                    {
                        if (x != maxCol) {
                            Qag(i, j) = Qag(i, j) + subQ[x](i, j);
                            nagg[i][j]++;
                        }
                    }
                }
            }

            for (int i = 0; i < col[maxCol]; i++)
            {
                for (int j = 0; j < col[maxCol]; j++)
                {
                    //cout << nagg[i][j] << " ,";
                    Qag(i, j) = Qag(i, j) / nagg[i][j];
                }
                //cout << endl;
            }

            cout << endl << "final AGG" << endl << Qag << endl;

            Ain = A.inverse();
            Bt = B.transpose();
            Q = Bt * Ain * B;

            MatrixXd Qago = Qorg.block(0, 0, col[maxCol], col[maxCol]);


            MatrixXd positiveDifferenceAgg = computePositiveDifference(Qago, Qag);
            cout << "Difference Matrix of Consolidated aggregate matrix with Positive Values:\n" << positiveDifferenceAgg << "\n";
            double avgPercentageErrorAgg = calculateAveragePercentageError(Qago, Qag);
            MatrixXd percentageErrorMatrixAgg = calculatePercentageErrorMatrix(Qago, Qag);
            cout << "Percentage Error Matrix of Consolidated aggregate matrix:\n" << percentageErrorMatrixAgg << "\n\n";
            cout << "Average Percentage Error of Consolidated aggregate matrix : " << avgPercentageErrorAgg << "%\n";
            double frobeniusNormError = calculateFrobeniusNormError(Qago, Qag);
            cout << "Frobenius Norm Error: " << frobeniusNormError << "\n";
            system("pause");
        }
       if(ch==3){
            system("CLS");
            int rank = 4;
            int size = 32;
            int rs = size - rank;
            //MatrixXd kernel(size, size); // Create a 16x16 MatrixXd

        // File handling to read matrix values from a file
            ifstream inputFile("matrix_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

            //inputFile >> rank >> size;


            MatrixXd kernel(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile >> kernel(i, j);
                }
            }

            inputFile.close();

            // Fill the matrix with some values (this is just an example)


            //MatrixXd kernel(size, size); // Create a 16x16 MatrixXd

        // File handling to read matrix values from a file
            ifstream inputFile1("matrix2_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile1.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

            // inputFile1 >> rank >> size;


            MatrixXd kernel1(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile1 >> kernel(i, j);
                }
            }

            inputFile1.close();

            //MatrixXd kernel(size, size); // Create a 16x16 MatrixXd

        // File handling to read matrix values from a file
            ifstream inputFile2("matrix3_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile2.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

            // inputFile2 >> rank >> size;


            MatrixXd kernel2(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile2 >> kernel(i, j);
                }
            }

            inputFile2.close();

            //MatrixXd kernel(size, size); // Create a 16x16 MatrixXd

        // File handling to read matrix values from a file
            ifstream inputFile3("matrix5_values.txt"); // Replace "matrix_values.txt" with your actual file name

            if (!inputFile3.is_open()) {
                cerr << "Error opening file." << endl;
                return 1;
            }

            //inputFile3 >> rank >> size;
            rs = size - rank;

            MatrixXd kernel3(size, size);

            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inputFile3 >> kernel(i, j);
                }
            }

            inputFile3.close();

            const int n = 4;
            array<MatrixXd, n> subB;
            array<MatrixXd, n> subA;
            array<MatrixXd, n> subAin;
            array<MatrixXd, n> subBt;
            array<MatrixXd, n> subQ;
            array<MatrixXd, n> subQor;
            array<MatrixXd, n> subQorg;



            // Define submatrices a, b, c, and d
            MatrixXd A = kernel.block(0, 0, rank, rank); // submatrix a (2x2)
            // submatrix b (2x14)

            subQorg[0] = kernel.block(rank, rank, rs, rs); // submatrix d (14x14)
            subQorg[1] = kernel1.block(rank, rank, rs, rs); // submatrix d (14x14)
            subQorg[2] = kernel2.block(rank, rank, rs, rs); // submatrix b (2x14)
            subQorg[3] = kernel3.block(rank, rank, rs, rs); //
            //MatrixXd Qorg = kernel.block(rank, rank, rs, rs); // submatrix d (14x14)




            int row[n];
            int col[n];
            for (int i = 0; i < n; i++)
            {
                cout << "enter dimensions of client " << i + 1 << endl;
                cin >> row[i];
                if (row[i] > rank)
                {
                    cout << "invalid dimensions";
                    return 0;
                }
                cin >> col[i];
                if (col[i] > rs)
                {
                    cout << "invalid dimensions";
                    return 0;
                }
            }
            subB[0] = kernel.block(0, rank, row[0], col[0]); // submatrix b (2x14)
            subB[1] = kernel1.block(0, rank, row[1], col[1]); // submatrix b (2x14)
            subB[2] = kernel2.block(0, rank, row[2], col[2]); // submatrix b (2x14)
            subB[3] = kernel3.block(0, rank, row[3], col[3]);

            /*subAin[0] = A.inverse();
            subBt[0] = subB[0].transpose();
            subQ[0] = subBt[0] * subAin[0] * subB[0];
            cout << "final of " << endl << subQ[0] << endl;
            */


            for (int i = 0; i < n; i++)
            {

                subA[i] = A.block(0, 0, row[i], row[i]);
                subAin[i] = subA[i].inverse();
                subBt[i] = subB[i].transpose();
                subQ[i] = subBt[i] * subAin[i] * subB[i];
                subQor[i] = subQorg[i].block(0, 0, col[i], col[i]);

                cout << "final of " << i + 1 << endl << subQ[i] << endl;
                MatrixXd positiveDifference = computePositiveDifference(subQor[i], subQ[i]);
                cout << "Difference Matrix of " << i + 1 << " client with Positive Values : \n" << positiveDifference << "\n";
                double avgPercentageError = calculateAveragePercentageError(subQor[i], subQ[i]);
                MatrixXd percentageErrorMatrix = calculatePercentageErrorMatrix(subQor[i], subQ[i]);
                cout << "Percentage Error Matrix of " << i + 1 << " client with :\n" << percentageErrorMatrix << "\n\n";
                cout << "Average Percentage Error of " << i + 1 << " client with :\n " << avgPercentageError << "%\n";
                double frobeniusNormError = calculateFrobeniusNormError(subQor[i], subQ[i]);
                cout << "Frobenius Norm Error: " << frobeniusNormError << "\n";

            }





            /*for (int i = 0; i < n; i++)
            {
                subBt[i] = subB[i].transpose();
                subQ[i] = subBt[i] * Ain * subB[i];
            }*/

            int maxCol = 0;
            for (int i = 1; i < n; ++i) {
                if (col[i] > col[maxCol]) {
                    maxCol = i;
                }
            }
            //cout << maxCol;
            MatrixXd Qag(col[maxCol], col[maxCol]);
            Qag = subQ[maxCol];

            int** nagg = new int* [col[maxCol]];
            for (int i = 0; i < col[maxCol]; ++i) {
                nagg[i] = new int[col[maxCol]];
            }

            for (int i = 0; i < col[maxCol]; ++i) {
                for (int j = 0; j < col[maxCol]; ++j) {
                    nagg[i][j] = 1;
                }
            }

            /*for (int i = 0; i < col[maxCol]; ++i) {
                for (int j = 0; j < col[maxCol]; ++j) {
                    cout<<nagg[i][j] <<", ";
                }
                cout << endl;
            }*/

            for (int x = 0; x < n; x++)
            {
                for (int i = 0; i < col[x]; i++)
                {
                    for (int j = 0; j < col[x]; j++)
                    {
                        if (x != maxCol) {
                            Qag(i, j) = Qag(i, j) + subQ[x](i, j);
                            nagg[i][j]++;
                        }
                    }
                }
            }

            cout << endl << "final AGG" << endl;
            for (int i = 0; i < col[maxCol]; i++)
            {
                for (int j = 0; j < col[maxCol]; j++)
                {
                    //cout << nagg[i][j] << " ,";
                    Qag(i, j) = Qag(i, j) / nagg[i][j];
                }
                //cout << endl;
            }
            MatrixXd QF = subQorg[0];
            for (int x = 1; x < n; x++)
            {
                for (int i = 0; i < col[maxCol]; i++)
                {
                    for (int j = 0; j < col[maxCol]; j++)
                    {
                        QF(i, j) = QF(i, j) + subQorg[x](i, j);

                    }
                }
            }

            cout << endl << "final AGG" << endl << Qag << endl;

            MatrixXd positiveDifference = computePositiveDifference(Qag, QF);
            cout << "Difference Matrix with Positive Values:\n" << positiveDifference << "\n";
            double avgPercentageError = calculateAveragePercentageError(Qag, QF);
            cout << "Average Percentage Error: " << avgPercentageError << "%\n";
            double frobeniusNormError = calculateFrobeniusNormError(Qag, QF);
            cout << "Frobenius Norm Error: " << frobeniusNormError << "\n";

           }
           if (ch == 4) {
               system("exit");
            }
           else {
               cout << "INVALID INPUT" << endl;
               exit;
           }

           goto label;
       /* case 4:
            system("CLS");
            cout << "\n\n\t\t\tThank You";
            break;
        default:cout << "\a";
        }
        cin.ignore();
        cin.get();
    } while (ch != 4);*/
    return 0;
}
