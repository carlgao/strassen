import java.util.Random;


public class Main 
{
	public static final int THRESHOLD = 2;

	public static void main (String[] args)
	{
		int DIMENSION = Integer.parseInt(args[0]);
		int d2 = DIMENSION + (DIMENSION & 1);
		
		int[][] a = randMatrix(d2, 0, 1);//{{1,0,0},{1,1,1},{0,0,1}};//
		int[][] b = randMatrix(d2, 0, 1);//{{0,0,1},{1,0,0},{1,1,1}};//randMatrix(2, 0, 1);

		int[][] c = new int[DIMENSION][DIMENSION];
		int[][] d = new int[DIMENSION][DIMENSION];
		if (1 == (DIMENSION & 1)) {
			int i = DIMENSION;
			for (int j = 0; j < d2; j++)
			{
				a[i][j] = 0;
				b[i][j] = 0;
			}
			int j = DIMENSION;
			for (i = 0; i < d2; i++)
			{
				a[i][j] = 0;
				b[i][j] = 0;
			}
			c = new int[d2][d2];
			d = new int[d2][d2];
		}
			
		

		
		mult2(a, b, c, DIMENSION);
		
		printMatrix(a);
		System.out.println("");
		printMatrix(b);
		System.out.println("");
		printMatrix(c);
		System.out.println("");

		strassen(a, b, d, DIMENSION);
		System.out.println("Strassen's: ");
		printMatrix(d);
		
		System.out.println("Actual: ");
		printMatrix(c);

		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				if (c[i][j] != d[i][j])
				{
					System.out.println("Different");
					return;
				}
			}
		}
		System.out.println("Same");
	}
	
	/*
	 * m1 = A B
	 *      C D
	 *      
	 * m2 = E F
	 *      G H 
	 * 
	 */
	public static void strassen(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		if (dim < THRESHOLD)
			mult(m1, m2, res, dim);
		else
		{
			System.out.println("lengths: " + m1.length + " " + m2.length);
			int dm = (dim+1)>>1; // account for padding
			int diff = dim - dm; //TODO: use dm everywhere instead of diff?

			int dm2 = dm + (dm & 1);

			int[][] temp1 = new int[dm2][dm2];
			int[][] temp2 = new int[dm2][dm2];


			// P1 = A(F-H)
			int[][] p1 = new int[dm2][dm2];		
			// copy a into temp1
			// copy f-h into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					System.out.println("P1: "+"m2["+(i)+","+(j+dm)+"]="+m2[i][j+dm]+"  m2["+(i+dm)+","+(j+dm)+"]="+m2[i+dm][j+dm]);
					temp1[i][j] = m1[i][j];
					temp2[i][j] = m2[i][j+dm] - m2[i+dm][j+dm];
				}
			}
			strassen(temp1, temp2, p1, dm);
			printMatrix(p1);

			// P2 = (A+B)H
			int[][] p2 = new int[dm2][dm2];		
			// copy a+b into temp1
			// copy h into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					System.out.println("P2: "+"m1["+(i)+","+(j)+"]="+m1[i][j]+"  m1["+(i)+","+(j+dm)+"]="+m1[i][j+dm]);
					temp1[i][j] = m1[i][j] + m1[i][j+dm];
					temp2[i][j] = m2[i+dm][j+dm];
				}
			}
			strassen(temp1, temp2, p2, dm);
			printMatrix(p2);


			// P3 = (C+D)E
			int[][] p3 = new int[dm2][dm2];		
			// copy c+d into temp1
			// copy e into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i+dm][j] + m1[i+dm][j+dm];
					temp2[i][j] = m2[i][j];
				}
			}
			strassen(temp1, temp2, p3, dm);
			printMatrix(p3);

			// P4 = D(G-E)
			int[][] p4 = new int[dm2][dm2];		
			// copy d into temp1
			// copy g-e into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i+dm][j+dm];
					temp2[i][j] = m2[i+dm][j] - m2[i][j];
				}
			}
			strassen(temp1, temp2, p4, dm);
			printMatrix(p4);

			// P5 = (A+D)(E+H)
			int[][] p5 = new int[dm2][dm2];
			// copy a+d into temp1, e+h into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j] + m1[i+dm][j+dm];
					temp2[i][j] = m2[i][j] + m2[i+dm][j+dm];
				}
			}		
			strassen(temp1, temp2, p5, dm);
			printMatrix(p5);

			// P6 = (B-D)(G+H)
			int[][] p6 = new int[dm2][dm2];
			// copy b-d into temp1, g+h into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j+dm] - m1[i+dm][j+dm];
					temp2[i][j] = m2[i+dm][j] + m2[i+dm][j+dm];
				}
			}		
			strassen(temp1, temp2, p6, dm);
			printMatrix(p6);

			// P7
			int[][] p7 = new int[dm2][dm2];
			// copy a-c into temp1, e+f into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j] - m1[i+dm][j];
					temp2[i][j] = m2[i][j] + m2[i][j+dm];
				}
			}		
			strassen(temp1, temp2, p7, dm);
			printMatrix(p7);

			// calculate result matrix
			System.out.println("calculating result matrix");
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					System.out.println("1: "+(p5[i][j] + p4[i][j] - p2[i][j] +p6[i][j]));
					res[i][j] = p5[i][j] + p4[i][j] - p2[i][j] +p6[i][j];
					
					System.out.println("2: "+(p1[i][j] + p2[i][j]));
					res[i][j+dm] = p1[i][j] + p2[i][j];
					
					System.out.println("3: "+(p3[i][j] + p4[i][j]));
					res[i+dm][j] = p3[i][j] + p4[i][j];
					
					System.out.println("4: "+(p5[i][j] + p1[i][j] - p3[i][j] -p7[i][j]));
					res[i+dm][j+dm] = p5[i][j] + p1[i][j] - p3[i][j] - p7[i][j];
				}
			}
			System.out.println("");
		}
	}
	
	public static int[][] randMatrix(int dim, int min, int max)
	{
		Random gen = new Random();
		
		int[][] mat = new int[dim][dim];
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				mat[i][j] = gen.nextInt(max-min+1)+min;
			}
		}
		return mat;
	}

	public static void printMatrix(int[][] m)
	{
		for (int i = 0; i < m.length; i++)
		{
			for (int j = 0; j < m[0].length; j++)
				System.out.print(m[i][j] + "  ");
			System.out.println("");
		}
	}
	
	// cache-efficient standard multiplication algorithm
	public static void mult2(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		for (int k = 0; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] = res[i][j] + (m1[i][k] * m2[k][j]);		
				}
	}
	
	public static void mult(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		System.out.println("dim: "+dim);
		for (int k = 0; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] = res[i][j] + (m1[i][k] * m2[k][j]);		
				}
	}
}

