#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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
                    tempX[i] += A[i*n+g] *X[g];
                }
            }
            tempX[i] /= -1* A[i*n+i];
        }
        norm = fabs(X[0]-tempX[0]);
        for(int h=0; h<n;h++)
        {
            if(fabs(X[h - tempX[h]]) > norm)
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

void MainWindow::on_action_Open_triggered()
{
    QString fileName = QFileDialog::getOpenFileName(this,tr("Open Linear System"),".",tr("Linear System Files (*.txt"));
    if(!fileName.isEmpty())
        loadFile(fileName);
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
    while(!in.atEnd())
    {
        QString line = in.readLine();
        QStringList lineList = line.split(" ");
        foreach(QString string, lineList)
        {
            bool ok=true;
            double value = string.toDouble(&ok);
            if(ok)
            {
                matrix.push_back(value);
                qDebug()<<value;
            }else qDebug() << "Error";
        }
    qDebug() << "end of line";
    }
    systemsFile.close();
}
