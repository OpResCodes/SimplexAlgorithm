using MathNet.Numerics.LinearAlgebra;
using NLog;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;



namespace SimpleSimplex
{
    public class Algorithm
    {

        private static Logger _log = LogManager.GetCurrentClassLogger();

        public void Optimize(LpModel model)
        {
            int m = model.A.RowCount;
            int n = model.A.ColumnCount;
            Iteration currentIteration = GetStartingBasis(model);
            while (true)
            {
                Matrix<double> B = GetColumnSubset(model, currentIteration.BasicColumns);
                Matrix<double> NonB = GetColumnSubset(model, currentIteration.NonBasicColumns);

                Matrix<double> B_inv = B.Inverse();
                Vector<double> Xb = B_inv.Multiply(model.b);

                _log.Info($"Iteration 0: Xb = {Xb.ToString()}");
                //reduced cost
                currentIteration.BasicCostVector = GetCostVector(model, currentIteration.BasicColumns);
                currentIteration.NonBasicCostVector = GetCostVector(model, currentIteration.NonBasicColumns);
                currentIteration.ObjectiveValue = currentIteration.BasicCostVector.ToRowMatrix().Multiply(Xb)[0];

                _log.Info($"Objective: {currentIteration.ObjectiveValue}");
                Matrix<double> cbInv = currentIteration.BasicCostVector.ToRowMatrix().Multiply(B_inv);
                Matrix<double> reducedCosts = cbInv.Multiply(NonB) - Matrix<double>.Build.DenseOfRowVectors(currentIteration.NonBasicCostVector);

                var min = reducedCosts.Row(0).Minimum();
                if (min >= 0)
                    break;
                var minId = reducedCosts.Row(0).MinimumIndex();//entering column etc, Taha 314
                int enteringColumn = currentIteration.NonBasicColumns[minId];
                _log.Info($"The Entering column is: {enteringColumn}");
                var bp = B_inv.Multiply(model.A.Column(enteringColumn));
                //check if unbounded.. if all entries negative or zero

                //ratio test to determine leaving vector
                //ratio test for positive denominator only -> set rest to NaN to ignore..
                bp.Map((a) =>
                {
                    if (a <= 0)
                        return double.NaN;
                    else
                        return a;
                }, bp);
                var feas = Xb.PointwiseDivide(bp);
                int leavingColumn = currentIteration.BasicColumns[feas.MinimumIndex()];
                _log.Info($"Leaving column: {leavingColumn}");
                //modify basis and continue
                currentIteration = ExchangeBasisColumns(model, currentIteration, enteringColumn, leavingColumn);
            }

            _log.Info($"Final Objective is: {Math.Round(currentIteration.ObjectiveValue, 2)}");

        }

        private Iteration ExchangeBasisColumns(LpModel model, Iteration previous, int entering, int leaving)
        {
            int[] nxtBasis = new int[previous.BasicColumns.Length];
            for (int i = 0; i < previous.BasicColumns.Length; i++)
            {
                if (previous.BasicColumns[i] == leaving)
                {
                    nxtBasis[i] = entering;
                }
                else
                {
                    nxtBasis[i] = previous.BasicColumns[i];
                }
            }
            int[] nxtNonBasis = new int[previous.NonBasicColumns.Length];
            for (int i = 0; i < previous.NonBasicColumns.Length; i++)
            {
                if (previous.NonBasicColumns[i] == entering)
                {
                    nxtNonBasis[i] = leaving;
                }
                else
                {
                    nxtNonBasis[i] = previous.NonBasicColumns[i];
                }
            }
            return new Iteration(nxtBasis, nxtNonBasis);
        }

        private Matrix<double> GetColumnSubset(LpModel model, int[] columns)
        {
            Vector<double>[] vecB = new Vector<double>[columns.Length];
            for (int i = 0; i < columns.Length; i++)
            {
                vecB[i] = model.A.Column(columns[i]);
            }
            return Matrix<double>.Build.DenseOfColumnVectors(vecB);
        }

        private Iteration GetStartingBasis(LpModel model)
        {
            int m = model.A.RowCount;
            int n = model.A.ColumnCount;
            int[] basicColumns = new int[m];
            int[] nonBasicColumns = new int[n - m];
            //first base
            for (int c = 0; c < n; c++)
            {
                if (c < (n - m))
                {
                    nonBasicColumns[c] = c;
                }
                else // c >= (n-m)
                {
                    basicColumns[c - (n - m)] = c;
                }
            }
            return new Iteration(basicColumns, nonBasicColumns);
        }

        Vector<double> GetCostVector(LpModel model, int[] columns)
        {
            double[] c = new double[columns.Length];
            for (int i = 0; i < columns.Length; i++)
            {
                c[i] = model.C[columns[i]];
            }
            return Vector<double>.Build.DenseOfArray(c);
        }

    }

    public struct Iteration
    {
        public Iteration(int[] basis, int[] nonBasis)
        {
            BasicColumns = basis;
            NonBasicColumns = nonBasis;
            BasicCostVector = null;
            NonBasicCostVector = null;
            ObjectiveValue = 0;
        }

        internal int[] BasicColumns { get; }

        internal int[] NonBasicColumns { get; }

        internal Vector<double> BasicCostVector { get; set; }

        internal Vector<double> NonBasicCostVector { get; set; }

        internal double ObjectiveValue { get; set; }

    }


    public class LpModel
    {
        public Vector<double> C { get; set; }

        public Matrix<double> A { get; set; }

        public Vector<double> b { get; set; }
    }
}
