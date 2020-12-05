package com.company;

import static java.lang.Math.*;

public class Solution {
    /**
     * Task solving procedures and boundary conditions
     */
    private final double xInit;
    private final double x0;
    private final int num;
    private final double h;
    private final double uInit;
    private final double uEnd;
    private final double error;
    private final double limit;
    private final Equation equation1;
    private final Equation equation2;

    public Solution(double xInit,
                    double xEnd,
                    double x0,
                    int num,
                    double uInit,
                    double uEnd,
                    double error,
                    int limit) {
        this.xInit = xInit;
        this.x0 = x0;
        this.num = num;
        this.uInit = uInit;
        this.uEnd = uEnd;
        this.error = error;
        this.limit = limit;
        h = (xEnd - xInit) / num;
        equation1 = new Equation(num);
        equation2 = new Equation(num);
    }

    // p-norm when p -> infinity
    private double norm() {
        int i = 0;
        double max = abs(equation1.get_u(i) - equation2.get_u(i));
        for (i = 1; i <= num; i++) {
            if (abs(equation1.get_u(i) - equation2.get_u(i)) > max) {
                max = abs(equation1.get_u(i) - equation2.get_u(i));
            }
        }
        return max;
    }

    // calculate solution of model task using TMA
    public void calc_mod() throws TmaException {
        int i = 0;
        equation1.init_grid(num, h);
        equation2.init_grid(num, h);
        equation2.double_grid();
        while (norm() / (pow(2,2) - 1) > error & i < limit | i == 0) {  // Runge rule
            equation1.tma_mod(xInit, x0);
            equation2.tma_mod(xInit, x0);
            equation1.double_grid();
            equation2.double_grid();
            i += 1;
        }
        if (i == limit)
            throw new TmaException("Accuracy is not achieved, increase the limit!");
    }

    // calculate solution of primary task using TMA
    public void calc() throws TmaException {
        int i = 0;
        equation1.init_grid(num, h);
        equation2.init_grid(num, h);
        equation2.double_grid();
        while (norm() / (pow(2,2) - 1) > error & i < limit | i == 0) {  // Runge rule
            equation1.tma(xInit, x0);
            equation2.tma(xInit, x0);
            equation1.double_grid();
            equation2.double_grid();
            i += 1;
        }
        if (i == limit)
            throw new TmaException("Accuracy is not achieved, increase the limit!");
    }

    // analytical solution of model task
    public void set_anal_mod() {
        double lambda1, lambda2;
        double mu1, mu2;
        double A11, A12, A21, A22;
        double B1, B2;
        double C1, C2, C3, C4;
        double x;
        double u;
        int la, lb;
        la = (int) floor((x0 - xInit) / h);  // position of the last dot before before discontinuity
        lb = la + 1;  // position of the first dot after discontinuity

        // lambda and mu when x < x0
        lambda1 = sqrt(Equation.q1(x0) / Equation.k1(x0));
        mu1 = Equation.f1(x0) / Equation.q1(x0);
        // lambda and mu when x > x0
        lambda2 = sqrt(Equation.q2(x0) / Equation.k2(x0));
        mu2 = Equation.f2(x0) / Equation.q2(x0);

        A11 = exp(-lambda1 * x0) - exp(lambda1 * x0);
        A12 = exp(lambda2 * (2 - x0)) - exp(lambda2 * x0);
        A21 = Equation.k1(x0) * lambda1 * (exp(lambda1 * x0) + exp(-lambda1 * x0));
        A22 = Equation.k2(x0) * lambda2 * (exp(lambda2 * (2 - x0)) + exp(lambda2 * x0));
        B1 = mu2 - mu1 + (mu1 - uInit) * exp(lambda1 * x0) - (mu2 - uEnd) * exp(lambda2 * (1 - x0));
        B2 = Equation.k1(x0) * lambda1 * (uInit - mu1) * exp(lambda1 * x0)
                + Equation.k2(x0) * lambda2 * (uEnd - mu2) * exp(lambda2 * (1 - x0));

        // constants when x < x0
        C1 = (((uInit - mu1) * A11 - B1) * A22 - ((uInit - mu1) * A21 - B2) * A12) / (A11 * A22 - A12 * A21);
        C2 = (B1 * A22 - B2 * A12) / (A11 * A22 - A12 * A21);
        // constants when x > x0
        C3 = (B2 * A11 - B1 * A21) / (A11 * A22 - A12 * A21);
        C4 = (uEnd - mu2) * exp(lambda2) - C3 * exp(2 * lambda2);

        for (int l = 1; l <= la; l++) {
            x = xInit + h * l;
            u = C1 * exp(lambda1 * x) + C2 * exp(-lambda1 * x) + mu1;
            equation2.set_u(l, u);
        }
        for (int l = num - 1; l >= lb; l--) {
            x = xInit + h * l;
            u = C3 * exp(lambda2 * x) + C4 * exp(-lambda2 * x) + mu2;
            equation2.set_u(l, u);
        }
    }

    public void init() {
        equation1.init_grid(num, h);
        equation1.set_u(0, uInit);
        equation1.set_u(num, uEnd);
        equation2.init_grid(num, h);
        equation2.set_u(0, uInit);
        equation2.set_u(num, uEnd);
    }

    public double get(int i) {
        return equation2.get_u(i);
    }

    public double get_prev(int i) {
        return equation1.get_u(i);
    }

    public double get_h() {
        return equation1.get_h();
    }
}
