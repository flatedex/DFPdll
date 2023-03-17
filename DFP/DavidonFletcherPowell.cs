using System;

namespace DFP
{
	class DavidonFletcherPowell
	{
		// функция поиска шага t
		private static double LineSearch(Func<double[], double> function, double[,] point, double[] direction, double epsilon)
		{
			double alpha = 1.0;
			double loBound = 0.0;
			double hiBound = Double.PositiveInfinity;
			double[] newPoint = FromMatrixToVector(point);
			while (true)
			{
				double[] candidate = Add(newPoint, Scale(alpha, direction));
				double functionValue = function(candidate);
				double loValue = function(Add(newPoint, Scale(loBound, direction)));
				if (functionValue > loValue + epsilon)
				{
					hiBound = alpha;
					alpha = (loBound + hiBound) / 2.0;
				}
				else
				{
					double hiValue = function(Add(newPoint, Scale(hiBound, direction)));
					if (functionValue > hiValue + epsilon)
					{
						loBound = alpha;
						alpha = (loBound + hiBound) / 2.0;
					}
					else
					{
						return alpha;
					}
				}
			}
		}

		// функция нахождения градиента
		private static double[] Gradient(Func<double[], double> function, double[,] point, double epsilon)
		{
			int n = point.Length;
			double[] newPoint = FromMatrixToVector(point);
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				double[] perturbation = new double[n];
				perturbation[i] = epsilon;
				result[i] = (function(Add(newPoint, perturbation)) - function(Subtract(newPoint, perturbation))) / (2 * epsilon);
			}
			return result;
		}

		// функция суммирования векторов
		private static double[] Add(double[] vector1, double[] vector2)
		{
			int n = vector1.Length;
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				result[i] = vector1[i] + vector2[i];
			}
			return result;
		}
		private static double[] Add(double[,] vector1, double[] vector2)
		{
			int n = vector1.Length;
			double[] result = new double[n];
			double[] newVector = FromMatrixToVector(vector1);
			for (int i = 0; i < n; i++)
			{
				result[i] = newVector[i] + vector2[i];
			}
			return result;
		}
		// функция вычитания векторов
		private static double[] Subtract(double[] vector1, double[] vector2)
		{
			int n = vector1.Length;
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				result[i] = vector1[i] - vector2[i];
			}
			return result;
		}
		private static double[] Subtract(double[] vector1, double[,] vector2)
		{
			int n = vector1.Length;
			double[] result = new double[n];
			double[] newVector = FromMatrixToVector(vector2);
			for (int i = 0; i < n; i++)
			{
				result[i] = vector1[i] - newVector[i];
			}
			return result;
		}

		// функция умножения вектора на скаляр
		private static double[] Scale(double scalar, double[] vector)
		{
			int n = vector.Length;
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				result[i] = scalar * vector[i];
			}
			return result;
		}

		//функция умножения матрицы на скаляр
		private static double[,] ScaleMatrix(double scalar, double[,] matrix)
		{
			int n = matrix.GetLength(0);
			int m = matrix.GetLength(1);
			double[,] result = new double[n, m];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					result[i, j] = matrix[i, j] * scalar;
				}
			}
			return result;
		}

		// функция умножения матрицы на вектор
		private static double[] MatrixOnVectorDotProduct(double[] vector, double[,] matrix)
		{
			int n = vector.Length;
			int m = matrix.GetLength(1);
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				double sum = 0.0;
				for (int j = 0; j < m; j++)
				{
					sum += vector[i] * matrix[i, j];
				}
				result[i] = sum;
			}
			return result;
		}

		// функция умножения матрицы на матрицу
		private static double[,] MatrixOnMatrixDotProduct(double[,] matrix1, double[,] matrix2)
		{
			int n = matrix1.GetLength(0);
			int m = matrix1.GetLength(1);
			int p = matrix2.GetLength(1);
			double[,] result = new double[n, p];
			for (int j = 0; j < p; j++)
			{
				for (int i = 0; i < n; i++)
				{
					double sum = 0.0;
					for (int k = 0; k < m; k++)
					{
						sum += matrix1[i, k] * matrix2[k, j];
					}
					result[i, j] = sum;
				}
			}
			return result;
		}
		// функция скалярного произведения векторов
		private static double VectorOnVectorDotProduct(double[] vector1, double[,] vector2)
		{
			int n = vector1.GetLength(0);
			int p = vector2.GetLength(1);
			double result = 0;
			for (int i = 0; i < p; i++)
			{
				result += vector1[i] * vector2[0, i];
			}
			return result;
		}
		// функция вычитания матриц
		private static double[,] MatrixSubtraction(double[,] matrix1, double[,] matrix2)
		{
			int n = matrix1.GetLength(0);
			int m = matrix1.GetLength(1);
			double[,] result = new double[n, m];
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < m; j++)
				{
					result[i, j] = matrix1[i, j] - matrix2[i, j];
				}
			}
			return result;
		}
		// функция превращения вектора [vector.Length] в матрицу [0, vector.Length]
		private static double[,] MatrixTranspose(double[] vector)
		{
			int n = vector.Length;
			double[,] result = new double[n, 1];
			for (int i = 0; i < n; i++)
			{
				result[i, 0] = vector[i];
			}
			return result;
		}
		// фукнция перевода вектора [vector.Length] в матрицу [vector.Length, 0]
		private static double[,] FromVectorToMatrix(double[] vector)
		{
			int n = vector.GetLength(0);
			double[,] result = new double[n, 1];
			for (int i = 0; i < n; i++)
			{
				result[i, 0] = vector[i];
			}
			return result;
		}
		private static double[] FromMatrixToVector(double[,] matrix)
		{
			int n = matrix.GetLength(0);
			double[] result = new double[n];
			for (int i = 0; i < n; i++)
			{
				result[i] = matrix[i, 0];
			}
			return result;
		}
		// функция создания единичной матрицы размерностью n
		private static double[,] IdentityMatrix(int n)
		{
			double[,] result = new double[n, n];
			for (int i = 0; i < n; i++)
			{
				result[i, i] = 1.0;
			}
			return result;
		}
		// функция для вычисления нормы вектора
		private static double Norm(double[] x)
		{
			double norm = 0;
			for (int i = 0; i < x.Length; i++)
			{
				norm += x[i] * x[i];
			}
			return Math.Sqrt(norm);
		}
		public static double[,] Minimize(Func<double[], double> f, double[,] x0, double accuracy)
		{
			// Начальное приближение
			double[,] x = x0;
			// Инициализируем единичную матрицу
			double[,] A0 = IdentityMatrix(x.Length);
			// Устанавливаем максимальное число итераций
			int maxIterations = 30000;
			// Устанавливаем пороговое значение для нормы градиента, т.н. точность
			double gradientTolerance = 1e-6;
			double newAccuracy = 1e-3;

			int k;
			for (k = 0; k < maxIterations; k++)
			{
				// Вычисляем градиент в точке x
				double[] grad = Gradient(f, x, gradientTolerance);

				// Проверяем, достигнута ли заданная точность
				if (Norm(grad) < newAccuracy)
				{
					System.Console.WriteLine(k.ToString());
					return x;
				}

				// Вычисляем направление спуска
				double[] d = MatrixOnVectorDotProduct(grad, A0);

				// Вычисляем длину шага по направлению градиента
				double t = LineSearch(f, x, d, accuracy);

				// Выполняем шаг оптимизации
				double[] xNew = Add(x, Scale(t, d));
				double[,] xNewMatrix = FromVectorToMatrix(xNew);
				// Вычисляем разность градиентов в точках xNew и x
				double[] gradNew = Gradient(f, xNewMatrix, gradientTolerance);
				double[] deltaGrad = Subtract(gradNew, grad);

				// Вычисляем разность аргументов в точках xNew и x
				double[] deltaX = Subtract(xNew, x);

				double[,] numerator1 = MatrixOnMatrixDotProduct(FromVectorToMatrix(deltaX), MatrixTranspose(deltaX));
				double[,] numerator2 = MatrixOnMatrixDotProduct(MatrixOnMatrixDotProduct(MatrixOnMatrixDotProduct(A0, FromVectorToMatrix(deltaGrad)), MatrixTranspose(deltaGrad)), A0);
				double denominetor1 = VectorOnVectorDotProduct(deltaGrad, MatrixTranspose(deltaX));
				double denominetor2 = VectorOnVectorDotProduct(deltaGrad, MatrixOnMatrixDotProduct(MatrixTranspose(deltaGrad), A0));

				// Вычисляем матрицу Ak+1
				double[,] A = MatrixSubtraction(ScaleMatrix(1 / denominetor1, numerator1), ScaleMatrix(1 / denominetor2, numerator2));

				// Обновляем x и H
				x = xNewMatrix;
				A0 = A;
			}
			System.Console.WriteLine("   " + k.ToString());
			// Если не удалось достичь заданной точности за максимальное число итераций,
			// возвращаем последнее

			return x;
		}
	}
}