package com.company;

import static java.lang.System.*;
import static java.lang.Math.*;

public class Main {
    public static void main(String[] args) {
        double xInit = 0;
        double xEnd = 1;
        double x0 = 1 / sqrt(2);
        int num = 10;
        double uInit = 0;
        double uEnd = 1;
        double error = 0.0001;
        int limit = 10;

        Solution solution = new Solution(
                xInit,
                xEnd,
                x0,
                num,
                uInit,
                uEnd,
                error,
                limit);
        Solution solution1 = new Solution(
                xInit,
                xEnd,
                x0,
                num,
                uInit,
                uEnd,
                error,
                limit);
        Solution solution2 = new Solution(
                xInit,
                xEnd,
                x0,
                num,
                uInit,
                uEnd,
                error,
                limit);

        solution.init();
        solution.set_anal_mod();
        solution1.init();
        solution1.calc_mod();
        out.println("Model task solution:");
        for (int i = 0; i <= num; i++) {
            out.println(solution1.get_prev(i) + " " + solution1.get(i) + " " + solution.get(i));
        }  // out: u_model(2h) u_model(h) u_anal
        out.println("Step: " + solution1.get_h());

        solution2.init();
        solution2.calc();
        out.println("Solution of the task with variable coefficients:");
        for (int i = 0; i <= num; i++) {
            out.println(solution2.get_prev(i) + " " + solution2.get(i));
        }  // out: u(2h) u(h)
        out.println("Step: " + solution2.get_h());
    }
}
