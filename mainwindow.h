#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QDebug>
#include <QFile>
#include <QFileDialog>
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

private:
    Ui::MainWindow *ui;
    QVector<double> matrix;
    QVector<double> rightRow;
    QVector<double> simpleIterations(QVector<double> A, QVector<double> B);
    QVector<double> successiveApproximations(QVector<double> A, QVector<double> B);
    QVector<double> zeydel(QVector<double> A, QVector<double> B);
    QVector<double> nekrasov(QVector<double> A, QVector<double> B);
    QVector<double> nevyazka(QVector<double>A,QVector<double>B,QVector<double>X);
    int maximumNumberOfIterations;
    bool loadFile(const QString &fileName);
    double eps;
};

#endif // MAINWINDOW_H
