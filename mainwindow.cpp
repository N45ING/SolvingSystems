#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    eps = 0.001;
    maximumNumberOfIterations = 50;
    stream.setString(&outputString);
    stream << qSetFieldWidth(5);
    stream << qSetRealNumberPrecision(4);
    QShortcut *fOpen = new QShortcut(QKeySequence("Ctrl+O"),this);
    QObject::connect(fOpen,SIGNAL(activated()),this,SLOT(on_action_Open_triggered()));
    QShortcut *Exit = new QShortcut(QKeySequence("Ctrl+E"),this);
    QObject::connect(Exit,SIGNAL(activated()),this,SLOT(on_action_Exit_triggered()));
}

MainWindow::~MainWindow()
{
    delete ui;
}
QVector<double> MainWindow::successiveApproximations(QVector<double> A, QVector<double> B)
{
    bool metka = false;
    int k = 0;
    int n = B.size();
    QVector<double> errorVector(B.size());
    QVector<double> X(B.size());
    QVector<double> Xk(B.size());
    for(int i=0;i<B.size(); i++)
    {
        X[i]=0.0;
        Xk[i]=0.0;
        errorVector[i]=0.0;
    }
    double g;
    while(metka == false && k <maximumNumberOfIterations)
    {
        for(int i=0; i<n; i++)
        {
            double s = 0;
            for(int j=0;j<n;j++)
            {
                s = s + A[i*n+j] * X[j];
            }
            Xk[i] = X[i] - s + B[i];
        }
        metka = true;
        for(int i=0;i<n;i++)
        {
            g=fabs(Xk[i]-X[i]);
            if(fabs(Xk[i])>1)
            {
                g = g/fabs(Xk[i]);
            }
            if(g>eps)
            {
                metka = false;
            }
        }
        k = k+1;
        for(int j=0;j<n;j++)
        {
            X[j]=Xk[j];
        }
    }
    if( k>= maximumNumberOfIterations)
    {
         qDebug() << "Error solving method" << endl; // change this message later
         errorMessageString = "Method is not convergent";
         return errorVector;
    }
    else
    {
        return X;
    }
}
QVector<double> MainWindow::simpleIterations(QVector<double> A, QVector<double> B)
{
    int n = B.size();
    QVector<double> errorVector(B.size());
    QVector<double> X(B.size());
    QVector<double> tempX(B.size());
    for(int i=0;i<B.size(); i++)
    {
        X[i]=0.0;
        tempX[i]=0.0;
        errorVector[i]=0.0;
    }
    double norm;
    int counter = 0;
    do
    {
        for(int i=0; i<n;i++)
        {
            tempX[i] = -1.0*B[i];
            for(int g=0;g<n;g++)
            {
                if(i!=g)
                {
                    tempX[i] += (A[i*n+g] *X[g]);
                }
            }
            tempX[i] /= -1* A[i*n+i];
        }
        norm = fabs(X[0]-tempX[0]);
        for(int h=0; h<n;h++)
        {
            if(fabs(X[h] - tempX[h]) > norm)
            {
                norm = fabs(X[h]-tempX[h]);
            }
            X[h] =tempX[h];
        }
        counter ++;
    }while(norm > eps && counter < maximumNumberOfIterations);
    if(counter >= maximumNumberOfIterations)
    {
        qDebug() << "Error solving method" << endl; // change this message later
        errorMessageString = "Method is not convergent";
        return errorVector;
    }
    else
    {
        return X;
    }
}
QVector<double> MainWindow::zeydel(QVector<double> A, QVector<double> B)
{
    int n = B.size();
    int counter =0;
    QVector<double> errorVector(B.size());
    QVector<double> X(B.size());
    QVector<double> Xk(B.size());
    for(int i=0;i<B.size(); i++)
    {
        X[i]=0.0;
        Xk[i]=0.0;
        errorVector[i]=0.0;
    }
    bool exitParam = true;
    do
    {
        for(int i=0; i<n;i++)
        {
            double var =0.0;
            for(int j=0;j<n;j++) // here might be a problem
            {
                if(j!=i) var+=(A[i*n+j]*X[j]);
            }
            Xk[i]=X[i];
            X[i]=(B[i]-var)/A[i*n+i];
        }
        exitParam = true;
        for(int i=0; i<n; i++)
        {
            if(fabs(X[i]-Xk[i]) >= eps)
            {
                exitParam = false;
            }
        }
        counter++;
    } while(!(exitParam)&&counter<maximumNumberOfIterations);
    if(counter == maximumNumberOfIterations )
    {
        qDebug() << "Error solving method" << endl; // change this message later
        errorMessageString = "Method is not convergent";
        return errorVector;
    }
    else
    {
        return X;
    }
}
QVector<double> MainWindow::nekrasov(QVector<double> A, QVector<double> B)
{
    int n = B.size();
    int counter =0;
    QVector<double> errorVector(B.size());
    QVector<double> X(B.size());
    QVector<double> Xk(B.size());
    for(int i=0;i<B.size(); i++)
    {
        X[i]=0.0;
        Xk[i]=0.0;
        errorVector[i]=0.0;
    }
    bool exitParam = true;
    do
    {
        for(int i=0; i<n; i++)
        {
            double s1 =0;
            for(int j=0;j<i;j++)
            {
                s1 = s1 + A[i*n+j] *X[j];
            }
            double s2 = 0;
            for(int l=i+1;l<n;l++)
            {
                s2 = s2+A[i*n+l] *Xk[l];
            }
            X[i]=1.0/A[i*n+i] * (B[i]-s1-s2);
        }
        exitParam = true;
        for(int i=0;i<n;i++)
        {
            if(fabs(X[i]-Xk[i]) >= eps)
            {
                exitParam = false;
            }
        }
        counter++;
        for(int i=0;i<n;i++)
        {
            Xk[i]=X[i];
        }
    }while(!(exitParam) && counter <maximumNumberOfIterations);
    if(counter >= maximumNumberOfIterations)
    {
        qDebug() << "Error solving method" << endl; // change this message later
        errorMessageString = "Method is not convergent";
        return errorVector;
    }
    else
    {
        return X;
    }
}
QVector<double> MainWindow::nevyazka(QVector<double> A, QVector<double> B, QVector<double> X)
{
    int n = B.size();
    QVector<double> nev(B.size());
    double s = 0;
    for(int i=0;i<n;i++)
    {
        s=0;
        for(int j=0;j<n;j++)
        {
            s= s+A[i*n+j] * X[j];
        }
        nev[i]=s-B[i];
    }
    return nev;
}
double MainWindow::nevyazkaSumm(QVector<double> A, QVector<double> B, QVector<double> X)
{
    int n = B.size();
    QVector<double> nev(B.size());
    double s = 0;
    for(int i=0;i<n;i++)
    {
        s=0;
        for(int j=0;j<n;j++)
        {
            s= s+A[i*n+j] * X[j];
        }
        nev[i]=s-B[i];
    }
    double nevyazkaSumm;
    for(int i=0;i <nev.size();i++)
    {
        nevyazkaSumm+=nev[i]*nev[i];
    }
    return sqrt(nevyazkaSumm);
}

void MainWindow::on_action_Open_triggered()
{
    ui->plainTextEdit->setPlainText("");
    ui->firstMethod->setEnabled(false);
    ui->secondMethod->setEnabled(false);
    ui->firstMethodOutput->setEnabled(false);
    ui->secondMethodOutput->setEnabled(false);
    ui->pushButton->setEnabled(false);
    ui->epsEdit->setEnabled(false);
    ui->itLimitEdit->setEnabled(false);
    ui->firstMethodOutput->setPlainText("");
    ui->secondMethodOutput->setPlainText("");
    stream.flush();
    outputString.clear();
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open Linear System"),".",tr("Linear System Files (*.txt"));
    if(!fileName.isEmpty())
        if(loadFile(fileName))
        {
            ui->firstMethod->setEnabled(true);
            ui->secondMethod->setEnabled(true);
            ui->firstMethodOutput->setEnabled(true);
            ui->pushButton->setEnabled(true);
            ui->epsEdit->setEnabled(true);
            ui->itLimitEdit->setEnabled(true);
        };
}
bool MainWindow::loadFile(const QString &fileName)
{
    QFile systemsFile(fileName);
    if(!systemsFile.open(QFile::ReadOnly | QFile::Text))
    {
        qDebug () << "could not open file for reading";
        return false;
    }
    QTextStream in(&systemsFile);
    QVector<QStringList> linesListVector;
    while(!in.atEnd())
    {
        QString line = in.readLine();
        QStringList lineList = line.split(" ");
        linesListVector.push_back(lineList);
    }
    bool canContinue=true;
    for(int i=1;i<linesListVector.size();i++)
    {
        if(linesListVector[i].size() != linesListVector[i-1].size()) canContinue = false;
    }
    if (linesListVector.size() != linesListVector[0].size()-1) canContinue = false;
    if(canContinue)
    {
        for(int i=0;i<linesListVector.size();i++)
        {
            for(int j=0;j<linesListVector[i].size();j++)
            {
                if(j==linesListVector[i].size()-1)
                {
                    bool ok = true;
                    QString string = linesListVector[i][j];
                    double value = string.toDouble(&ok);
                    if(ok)
                    {
                        rightRow.push_back(value);
                        //qDebug()<<value;
                        stream <<" "<< value << endl;
                    }else
                    {
                        qDebug() << "Error";
                        stream.flush();
                        stream << "Error input, please, check txt file." ;
                        printf(outputString);
                        systemsFile.close();
                        return false;
                    }
                }else
                {
                    bool ok = true;
                    QString string = linesListVector[i][j];
                    double value = string.toDouble(&ok);
                    if(ok)
                    {
                        matrix.push_back(value);
                        //qDebug()<<value;
                        stream << value<<"  ";
                    }else
                    {
                        qDebug() << "Error";
                        stream.flush();
                        stream << "Error input" ;
                        printf(outputString);
                        systemsFile.close();
                        return false;
                    }
                 }
            }
        }
        printf(outputString);
    }else
    {
        qDebug() << "error input";
        stream.flush();
        stream << "Error input" ;
        printf(outputString);
        systemsFile.close();
        return false;
    }
    /*qDebug() << simpleIterations(matrix,rightRow);
    qDebug() << nevyazka(matrix,rightRow,simpleIterations(matrix,rightRow));*/
    systemsFile.close();
    return true;
}
void MainWindow::printf(QString string)
{
    ui->plainTextEdit->setPlainText(string);
}
void MainWindow::on_action_Exit_triggered()
{
    MainWindow::close();
}

void MainWindow::on_epsEdit_valueChanged(double arg1)
{
    eps = arg1;
}

void MainWindow::on_itLimitEdit_valueChanged(int arg1)
{
    maximumNumberOfIterations = arg1;
}

void MainWindow::on_secondMethod_currentIndexChanged(const QString &arg1)
{
    if (ui->secondMethod->currentIndex()!=0)
    {
        ui->secondMethodOutput->setEnabled(true);
    }
    else
    {
        ui->secondMethodOutput->setEnabled(false);
    }
}

void MainWindow::on_pushButton_clicked()
{
    ui->firstMethodOutput->setPlainText("");
    ui->secondMethodOutput->setPlainText("");
    errorMessageString="";
    QVector<double> X1;
    QVector<double> nev1;
    double nevSumm1;
    double nevSumm2;
    QVector<double> X2;
    QVector<double> nev2;
    if(ui->firstMethod->currentText()=="Successive Approximation")
    {
        X1=successiveApproximations(matrix,rightRow);
        nev1=nevyazka(matrix,rightRow,X1);
        nevSumm1=nevyazkaSumm(matrix,rightRow,X1);
    }
    if(ui->firstMethod->currentText()=="Simple Iteration")
    {
        X1=simpleIterations(matrix,rightRow);
        nev1=nevyazka(matrix,rightRow,X1);
        nevSumm1=nevyazkaSumm(matrix,rightRow,X1);
    }
    if(ui->firstMethod->currentText()=="Zeidel's Method")
    {
        X1=zeydel(matrix,rightRow);
        nev1=nevyazka(matrix,rightRow,X1);
        nevSumm1=nevyazkaSumm(matrix,rightRow,X1);
    }
    if(ui->firstMethod->currentText()=="Nekrasov's Method")
    {
        X1=nekrasov(matrix,rightRow);
        nev1=nevyazka(matrix,rightRow,X1);
        nevSumm1=nevyazkaSumm(matrix,rightRow,X1);
    }

    if(ui->secondMethod->currentText()=="Successive Approximation")
    {
        X2=successiveApproximations(matrix,rightRow);
        nev2=nevyazka(matrix,rightRow,X2);
        nevSumm2=nevyazkaSumm(matrix,rightRow,X2);
    }
    if(ui->secondMethod->currentText()=="Simple Iteration")
    {
        X2=simpleIterations(matrix,rightRow);
        nev2=nevyazka(matrix,rightRow,X2);
        nevSumm2=nevyazkaSumm(matrix,rightRow,X2);
    }
    if(ui->secondMethod->currentText()=="Zeidel's Method")
    {
        X2=zeydel(matrix,rightRow);
        nev2=nevyazka(matrix,rightRow,X2);
        nevSumm2=nevyazkaSumm(matrix,rightRow,X2);
    }
    if(ui->secondMethod->currentText()=="Nekrasov's Method")
    {
        X2=nekrasov(matrix,rightRow);
        nev2=nevyazka(matrix,rightRow,X2);
        nevSumm2=nevyazkaSumm(matrix,rightRow,X2);
    }
    QTextStream firstStream;
    QTextStream secondStream;
    QString firstString;
    QString secondString;
    firstStream.setString(&firstString);
    secondStream.setString(&secondString);
    if(errorMessageString.isEmpty())
    {
        firstStream << "Solution Vector:" << endl;
        for(int i=0;i<X1.size();i++)
        {
            firstStream << qSetRealNumberPrecision(5)<< X1[i] << " ";
        }
        firstStream << endl;
        firstStream << "Mistake Vector:" << endl;
        for(int i=0;i<nev1.size();i++)
        {
            firstStream << qSetRealNumberPrecision(5) << nev1[i] << " ";
        }
        firstStream << endl;
        firstStream << "Residual: ";
        firstStream << nevSumm1;
        firstStream << endl;
        ui->firstMethodOutput->setPlainText(firstString);
    } else
    {
        ui->firstMethodOutput->setPlainText(errorMessageString);
    }
    if(ui->secondMethod->currentText()!="Don't Compare" && errorMessageString.isEmpty())
    {
        secondStream << "Solution Vector:" << endl;
        for(int i=0;i<X2.size();i++)
        {
            secondStream << qSetRealNumberPrecision(5)<< X2[i] << " ";
        }
        secondStream<<endl;
        secondStream << "Mistake Vector:" << endl;
        for(int i=0;i<nev2.size();i++)
        {
            secondStream << qSetRealNumberPrecision(5) << nev2[i] << " ";
        }
        secondStream << endl;
        secondStream << "Residual: ";
        secondStream << nevSumm2;
        secondStream << endl;
        ui->secondMethodOutput->setPlainText(secondString);
        if(nevSumm2<nevSumm1)
        {
            QMessageBox::information(this,"Accuracy","Second method is more accurate");
        }else if(nevSumm2==nevSumm1)
        {
            QMessageBox::information(this,"Accuracy","Same Accuracy on both methods");
        }else
        {
            QMessageBox::information(this,"Accuracy","First method is more accurate");
        }
    }else if (ui->secondMethod->currentText()!="Don't Compare")
    {
        ui->secondMethodOutput->setPlainText(errorMessageString);
    }
}

