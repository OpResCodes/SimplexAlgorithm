using MathNet.Numerics.LinearAlgebra;
using NLog;
using System;

namespace SimplexTester
{
    class Program
    {
        static void Main(string[] args)
        {
            StartLogger();
            SimpleSimplex.LpModel lpModel = new SimpleSimplex.LpModel();
            lpModel.A = Matrix<double>.Build.DenseOfArray(new double[,] { 
                {  6,4,1,0,0,0},
                {  1,2,0,1,0,0},
                { -1,1,0,0,1,0},
                {  0,1,0,0,0,1}
            });
            lpModel.C = Vector<double>.Build.DenseOfArray(new double[] { 5,4,0,0,0,0});
            lpModel.b = Vector<double>.Build.DenseOfArray(new double[] { 24,6,1,2 });


            SimpleSimplex.Algorithm algorithm = new SimpleSimplex.Algorithm();
            algorithm.Optimize(lpModel);
            Console.WriteLine("key");
            Console.ReadKey();

        }

        private static void StartLogger()
        {
            //logging
            var config = new NLog.Config.LoggingConfiguration();
            var traceTarget = new NLog.Targets.TraceTarget();
            var consoleTarget = new NLog.Targets.ConsoleTarget();
            config.AddRule(LogLevel.Trace, LogLevel.Fatal, traceTarget);
            config.AddRule(LogLevel.Trace, LogLevel.Fatal, consoleTarget);
            LogManager.Configuration = config;
            LogManager.GetCurrentClassLogger().Info("Application started.");
        }
    }
}
