import java.util.Random;


public class Main 
{
	public static final int THRESHOLD = 2;
	public static final int DIMENSION = 16;

	public static void main (String[] args)
	{
		int[][] a = randMatrix(DIMENSION, 0, 1);//{{1,0,0},{1,1,1},{0,0,1}};//
		int[][] b = randMatrix(DIMENSION, 0, 1);//{{0,0,1},{1,0,0},{1,1,1}};//randMatrix(2, 0, 1);

		int[][] c = new int[DIMENSION][DIMENSION];
		int[][] d = new int[DIMENSION][DIMENSION];

		mult2(a, b, c, DIMENSION);
		
		printMatrix(a);
		System.out.println("");
		printMatrix(b);
		System.out.println("");
		printMatrix(c);
		System.out.println("");

		strassen(a, b, d, 0, 0, 0, 0, DIMENSION);
		System.out.println("Strassen's: ");
		printMatrix(d);
		
		System.out.println("Actual: ");
		printMatrix(c);

	}
	
	/*
	 * m1 = A B
	 *      C D
	 *      
	 * m2 = E F
	 *      G H 
	 * 
	 */
	public static void strassen(int[][] m1, int[][] m2, int[][] res, int or1, int oc1, int or2, int oc2, int dim)
	{
		if (dim < THRESHOLD)
			mult(m1, m2, res, or1, oc1, or2, oc2, dim);
		else
		{
			int dm = (dim+1)>>1; // account for padding
			int diff = dim - dm;


			int[][] temp1 = new int[dm][dm];
			int[][] temp2 = new int[dm][dm];


			// P1 = A(F-H)
			int[][] p1 = new int[dm][dm];		
			// copy f-h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					System.out.println("P1: "+"m2["+(i+or2)+","+(j+dm+oc2)+"]="+m2[i+or2][j+dm+oc2]+"  m2["+(i+dm+or2)+","+(j+dm+oc2)+"]="+m2[i+dm+or2][j+dm+oc2]);
					temp2[i][j] = m2[i+or2][j+dm+oc2] - m2[i+dm+or2][j+dm+oc2];
				}
			}
			strassen(m1, temp2, p1, or1, oc1, 0, 0, dm);
			printMatrix(p1);

			// P2 = (A+B)H
			int[][] p2 = new int[dm][dm];		
			// copy a+b into temp1
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					System.out.println("P2: "+"m1["+(i+or1)+","+(j+oc1)+"]="+m1[i+or1][j+oc1]+"  m1["+(i+or1)+","+(j+dm+oc1)+"]="+m1[i+or1][j+dm+oc1]);
					temp1[i][j] = m1[i+or1][j+oc1] + m1[i+or1][j+dm+oc1];
				}
			}
			strassen(temp1, m2, p2, 0, 0, or2 + dm, oc2 + dm, dm);
			printMatrix(p2);


			// P3 = (C+D)E
			int[][] p3 = new int[dm][dm];		
			// copy c+d into temp1
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i+dm+or1][j+oc1] + m1[i+dm+or1][j+dm+oc1];
				}
			}
			strassen(temp1, m2, p3, 0, 0, or2, oc2, dm);
			printMatrix(p3);

			// P4 = D(G-E)
			int[][] p4 = new int[dm][dm];		
			// copy g-e into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp2[i][j] = m2[i+dm+or2][j+oc2] - m2[i+or2][j+oc2];
				}
			}
			strassen(m1, temp2, p4, or1+dm, oc1+dm, 0, 0, dm);
			printMatrix(p4);

			// P5 = (A+D)(E+H)
			int[][] p5 = new int[dm][dm];
			// copy a+d into temp1, e+h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i+or1][j+oc1] + m1[i+dm+or1][j+dm+oc1];
					temp2[i][j] = m2[i+or2][j+oc2] + m2[i+dm+or2][j+dm+oc2];
				}
			}		
			strassen(temp1, temp2, p5, 0, 0, 0, 0, dm);
			printMatrix(p5);

			// P6 = (B-D)(G+H)
			int[][] p6 = new int[dm][dm];
			// copy b-d into temp1, g+h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i+or1][j+dm+oc1] - m1[i+dm+or1][j+dm+oc1];
					temp2[i][j] = m2[i+dm+or2][j+oc2] + m2[i+dm+or2][j+dm+oc2];
				}
			}		
			strassen(temp1, temp2, p6, 0, 0, 0, 0, dm);
			printMatrix(p6);

			// P7
			int[][] p7 = new int[dm][dm];
			// copy a-c into temp1, e+f into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i+or1][j+oc1] - m1[i+dm+or1][j+oc1];
					temp2[i][j] = m2[i+or2][j+oc2] + m2[i+or2][j+dm+oc2];
				}
			}		
			strassen(temp1, temp2, p7, 0, 0, 0, 0, dm);
			printMatrix(p7);

			// calculate result matrix
			System.out.println("calculating result matrix");
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
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
	
	public static void mult(int[][] m1, int[][] m2, int[][] res, int or1, int oc1, int or2, int oc2, int dim)
	{
		System.out.println("dim: "+dim);
		System.out.println(or1+" "+oc1+" "+or2+" "+oc2+" ");
		for (int k = 0; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] = res[i][j] + (m1[i+or1][k+oc1] * m2[k+or2][j+oc2]);		
				}
	}
}

