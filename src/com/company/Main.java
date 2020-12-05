package com.company;

import static java.lang.System.*;
import static java.lang.Math.*;
import java.io.*;

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

        File file = new File("Output.txt");
        try (FileWriter writer = new FileWriter(file, false)) {
            try {
                writer.write("Model task solution:\n");
                solution.init();
                solution.set_anal_mod();
                solution1.init();
                solution1.calc_mod();
                for (int i = 0; i <= num; i++) {
                    writer.write(solution1.get_prev(i) + " " + solution1.get(i) + " " + solution.get(i) + "\n");
                }  // out: u_model(2h) u_model(h) u_anal
                writer.write("Step: " + solution1.get_h() + "\n");

                writer.write("Solution of the task with variable coefficients:\n");
                solution2.init();
                solution2.calc();
                for (int i = 0; i <= num; i++) {
                    writer.write(solution2.get_prev(i) + " " + solution2.get(i) + "\n");
                }  // out: u(2h) u(h)
                writer.write("Step: " + solution2.get_h() + "\n");
            }
            catch (TmaException ex) {
                writer.write(ex.getMessage());
            }
        }
        catch (IOException ex) {
            out.println(ex.getMessage());
        }
    }
}
