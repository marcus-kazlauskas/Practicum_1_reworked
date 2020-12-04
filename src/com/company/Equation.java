package com.company;

import static java.lang.Math.*;

public class Equation {
    /**
     * Solution grid function, task condition functions and solution algorithm (TMA)
     */
    private int num;  // number of line segments
    private double h;  // step of calculating
    private final double[] u;  // grid solution of the equation

    public Equation(int num) {
        u = new double[num + 1];
    }

    // setter and getter
    public void set_u(int i, double u) {
        this.u[i] = u;
    }
    public double get_u(int i) {
        return u[i];
    }

    // grid
    public void init_grid(int num, double h) {
        this.num = num;
        this.h = h;
    }
    public void double_grid() {
        num *= 2;
        h /= 2;
    }
    public double get_h() {
        return h;
    }

    // k(x), q(x), f(x) when x < x0
    public static double k1(double x) {
        return pow(x, 2) + 0.5;
    }
    public static double q1(double x) {
        return exp(-pow(x, 2));
    }
    public static double f1(double x) {
        return cos(x);
    }

    // k(x), q(x), f(x) when x > x0
    public static double k2(double x) {
        return pow(x, 2) + 0.5;
    }
    public static double q2(double x) {
        return 1;
    }
    public static double f2(double x) {
        return 1;
    }

    // coefficient a(x), b(x), c(x), d(x) in TMA when x < x0
    // alpha(x), beta(x) are calculated in tma() from xInit to x0 with step h
    private double a1(double x) {
        return k1(x + h / 2);
    }
    private double b1(double x) {
        return -(k1(x + h / 2) + k1(x - h / 2) + q1(x) * pow(h, 2));
    }
    private double c1(double x) {
        return k1(x - h / 2);
    }
    private double d1(double x) {
        return -f1(x) * pow(h, 2);
    }

    // coefficients a(x), b(x), c(x), d(x) in TMA when x > x0
    // alpha(x), beta(x) are calculated in tma() from xEnd to x0 with step h
    private double a2(double x) {
        return k2(x + h / 2);
    }
    private double b2(double x) {
        return -(k2(x + h / 2) + k2(x - h / 2) + q2(x) * pow(h, 2));
    }
    private double c2(double x) {
        return k2(x - h / 2);
    }
    private double d2(double x) {
        return -f2(x) * pow(h, 2);
    }

    // a(x0), b(x0), c(x0), d(x0) in TMA for model task when x < x0
    // alpha(x), beta(x) are calculated in tma_mod() from xInit to x0 with step h
    private double a1_mod(double x0) {
        return k1(x0);
    }
    private double b1_mod(double x0) {
        return -2 * k1(x0) - q1(x0) * pow(h, 2);
    }
    private double c1_mod(double x0) {
        return k1(x0);
    }
    private double d1_mod(double x0) {
        return -f1(x0) * pow(h, 2);
    }

    // a(x0), b(x0), c(x0), d(x0) in TMA for model task when x > x0
    // alpha(x), beta(x) are calculated in tma_mod() from xEnd to x0 with step h
    private double a2_mod(double x0) {
        return k2(x0);
    }
    private double b2_mod(double x0) {
        return -2 * k2(x0) - q2(x0) * pow(h, 2);
    }
    private double c2_mod(double x0) {
        return k2(x0);
    }
    private double d2_mod(double x0) {
        return -f2(x0) * pow(h, 2);
    }

    // TMA for model task
    public void tma_mod(double xInit, double x0) {
        double[] alpha = new double[num + 1];
        double[] beta = new double[num + 1];
        double uBuf;
        int la, lb;
        int n = num / (u.length - 1);

        alpha[0] = 0;
        beta[0] = u[0];
        alpha[num] = 0;
        beta[num] = u[u.length - 1];
        la = (int) floor((x0 - xInit) / h);                   // точка la до точки разрыва
        lb = la + 1;                                  // точка lb с другой стороны от разрыва

        for (int l = 1; l < la; l++) {                                                   // прогонка от 0 до la
            alpha[l] = -a1_mod(x0) / (b1_mod(x0) + c1_mod(x0) * alpha[l - 1]);
            beta[l] = (d1_mod(x0) - c1_mod(x0) * beta[l - 1]) / (b1_mod(x0) + c1_mod(x0) * alpha[l - 1]);
        }
        for (int l = num - 1; l > lb; l--) {                                              // прогонка от 1 до lb
            alpha[l] = -c2_mod(x0) / (b2_mod(x0) + a2_mod(x0) * alpha[l + 1]);
            beta[l] = (d2_mod(x0) - a2_mod(x0) * beta[l + 1]) / (b2_mod(x0) + a2_mod(x0) * alpha[l + 1]);
        }

        uBuf = (k1(x0) * beta[la - 1] + k2(x0) * beta[lb + 1])
                / (k1(x0) * (1 - alpha[la - 1]) + k2(x0) * (1 - alpha[lb + 1]));   // u(la) = u(lb)
        for (int l = la - 1; l >= 1; l--) {                                                  // обратная прогонка от u(la) до u(0+0)
            uBuf = alpha[l] * uBuf + beta[l];
            if (l % n == 0) {                                                          // запись одной из 11-ти точек
                u[l / n] = uBuf;
            }
        }

        uBuf = (k1(x0) * beta[la - 1] + k2(x0) * beta[lb + 1])
                / (k1(x0) * (1 - alpha[la - 1]) + k2(x0) * (1 - alpha[lb + 1]));   // u(lb) = u(la)
        for (int l = lb + 1; l <= num - 1; l++) {                                             // обратная прогонка от u(lb) до (1-0)
            uBuf = alpha[l] * uBuf + beta[l];
            if (l % n == 0) {                                                          // запись одной из 11-ти точек
                u[l / n] = uBuf;
            }
        }
    }

    // TMA for primary task
    public void tma(double xInit, double x0) {
        double[] alpha = new double[num + 1];
        double[] beta = new double[num + 1];
        double x;
        double uBuf;
        int la, lb;
        double xa, xb;
        int n = num / (u.length - 1);

        alpha[0] = 0;
        beta[0] = u[0];
        alpha[num] = 0;
        beta[num] = u[u.length - 1];
        la = (int) floor((x0 - xInit) / h);                   // точка la до точки разрыва
        lb = la + 1;                                  // точка lb с другой стороны от разрыва
        xa = xInit + h * la;
        xb = xInit + h * lb;

        for (int l = 1; l < la; l++) {                                                   // прогонка от 0 до la
            x = xInit + h * l;
            alpha[l] = -a1(x) / (b1(x) + c1(x) * alpha[l - 1]);
            beta[l] = (d1(x) - c1(x) * beta[l - 1]) / (b1(x) + c1(x) * alpha[l - 1]);
        }
        for (int l = num - 1; l > lb; l--) {                                              // прогонка от 1 до lb
            x = xInit + h * l;
            alpha[l] = -c2(x) / (b2(x) + a2(x) * alpha[l + 1]);
            beta[l] = (d2(x) - a2(x) * beta[l + 1]) / (b2(x) + a2(x) * alpha[l + 1]);
        }

        uBuf = (k1(xa) * beta[la - 1] + k2(xb) * beta[lb + 1])
                / (k1(xa) * (1 - alpha[la - 1]) + k2(xb) * (1 - alpha[lb + 1]));   // u(la) = u(lb)
        for (int l = la - 1; l >= 1; l--) {                                                  // обратная прогонка от u(la) до u(0+0)
            uBuf = alpha[l] * uBuf + beta[l];
            if (l % n == 0) {                                                          // запись одной из 11-ти точек
                u[l / n] = uBuf;
            }
        }

        uBuf = (k1(xa) * beta[la - 1] + k2(xb) * beta[lb + 1])
                / (k1(xa) * (1 - alpha[la - 1]) + k2(xb) * (1 - alpha[lb + 1]));   // u(lb) = u(la)
        for (int l = lb + 1; l <= num - 1; l++) {                                             // обратная прогонка от u(lb) до (1-0)
            uBuf = alpha[l] * uBuf + beta[l];
            if (l % n == 0) {                                                          // запись одной из 11-ти точек
                u[l / n] = uBuf;
            }
        }
    }
}
