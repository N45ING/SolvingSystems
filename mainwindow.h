#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
#include <QShortcut>
#include <QDrag>
#include <QMimeData>
#include <QMessageBox>
#include <math.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    void on_action_Open_triggered();

    void on_action_Exit_triggered();

    void on_epsEdit_valueChanged(double arg1);

    void on_itLimitEdit_valueChanged(int arg1);

    void on_secondMethod_currentIndexChanged(const QString &arg1);

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
    QVector<double> matrix;
    QVector<double> rightRow;
    QVector<double> simpleIterations(QVector<double> A, QVector<double> B);
    QVector<double> successiveApproximations(QVector<double> A, QVector<double> B);
    QVector<double> zeydel(QVector<double> A, QVector<double> B);
    QVector<double> nekrasov(QVector<double> A, QVector<double> B);
    double nevyazkaSumm(QVector<double>A,QVector<double>B,QVector<double>X);
    QVector<double> nevyazka(QVector<double> A, QVector<double> B, QVector<double> X);
    int maximumNumberOfIterations;
    bool loadFile(const QString &fileName);
    double eps;
    QTextStream stream;
    QString outputString;
    QString errorMessageString;
    void printf(QString string);
};

#endif // MAINWINDOW_H
