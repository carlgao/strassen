import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;


public class TryAgain 
{
	static int DIMENSION;
	static int THRESHOLD;
	static int TRIALS;

	static HashMap<Integer, int[][]> P;
	static HashMap<Integer, int[][]> T1;
	static HashMap<Integer, int[][]> T2;

	public static void main(String[] args)
	{
		DIMENSION = 120; //Integer.parseInt(args[1]);
		THRESHOLD = DIMENSION; // Integer.parseInt(args[2]);
		TRIALS = 1000;
		
		int[][] a = randMatrix(DIMENSION, 0, 13341);
		int[][] b = randMatrix(DIMENSION, 0, 4);
		// readMatrices(args[3], a, b);

		initStrasMem();

		timeMult(a, b);
		System.out.println("");
		timeStrassen(a, b);
		
		System.out.println("");
		/*
		if (equals(c,d, DIMENSION))
			System.out.println("YAYYYY");
		else
			System.out.println("FML");*/

	}
	
	public static void timeStrassen(int[][] a, int[][] b)
	{
		long t = 0;
		long ttemp = -1;
		
		// discard first 10 results because the program seems to be slow when "warming up" initially, 
		// probably due to caching effects
		for (int i = 0; i < 10; i++)
		{
			int[][] c = new int[DIMENSION][DIMENSION];
			strassen(a, b, c, DIMENSION);
		}

		// average over n=TRIALS matrix multiplications
		for (int i = 0; i < TRIALS; i++)
		{
			int[][] c = new int[DIMENSION][DIMENSION];
			ttemp = System.nanoTime();
			strassen(a, b, c, DIMENSION);
			ttemp = System.nanoTime() - ttemp;
			t += ttemp;
		}
		t /= TRIALS;

		System.out.println("STRASSEN");
		System.out.println("Threshold: "+THRESHOLD);
		System.out.println("Dimension: "+DIMENSION);
		System.out.println("Average computation time: "+ t + " ns");
	}
	
	public static void timeMult(int[][] a, int[][] b)
	{
		long t = 0;
		long ttemp = -1;
		
		// discard first 10 results because the program seems to be slow when "warming up" initially, 
		// probably due to caching effects
		for (int i = 0; i < 10; i++)
		{
			int[][] c = new int[DIMENSION][DIMENSION];
			mult(a, b, c, DIMENSION);
		}
		// average over n=TRIALS matrix multiplications
		for (int i = 0; i < TRIALS; i++)
		{
			int[][] c = new int[DIMENSION][DIMENSION];

			ttemp = System.nanoTime();
			mult(a, b, c, DIMENSION);
			ttemp = System.nanoTime() - ttemp;
			t += ttemp;
		}
		t /= TRIALS;

		System.out.println("CONVENTIONAL MULTIPLICATION");
		System.out.println("Dimension: "+DIMENSION);
		System.out.println("Average computation time: "+ t + " ns");
	}
	
	/*
	 * Initialize storage for all dimensions of Strassen's recursion.
	 */
	public static void initStrasMem()
	{
		P = new HashMap<Integer, int[][]>();
		T1 = new HashMap<Integer, int[][]>();
		T2 = new HashMap<Integer, int[][]>();

		int size = DIMENSION;

		while (size >= THRESHOLD)
		{
			size = (size + 1) >> 1;
			//System.out.println("adding "+ size);
			P.put(size, new int[size][size]);
			T1.put(size, new int[size][size]);
			T2.put(size, new int[size][size]);
		}
	}
	
	/*
	 * m1 = A B             m2 = E F         res = (AE+BG)  (AF+BH)
	 *      C D                  G H               (CE+DG)  (CF+DH)
	 *      
	 * Combinations:
	 *    AE+BG = P5+P4-P2+P6
	 *    AF+BH = P1+P2
	 *    CE+DG = P3+P4
	 *    CF+DH = P5+P1-P3-P7
	 *    
	 */
	public static void strassen(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		if (dim < THRESHOLD)
			mult(m1, m2, res, dim);
		
		else
		{
			int dhalf = (dim+1) >> 1;
			boolean odd = (dim & 1) == 1;
		
			int[][] t1 = T1.get(dhalf);
			if (t1 == null)
			{
				System.out.println(dhalf+ " fail");
				System.exit(0);
			}
			int[][] t2 = T2.get(dhalf);
			int[][] p = P.get(dhalf);
			
			
			// P1 = A(F-H) ===========================
			// copy A into t1
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t1[i][j] = m1[i][j];
			// copy F into temp2
			for (int i = 0; i < dhalf; i++)
			{
				for (int j = dhalf; j < dim; j++)
					t2[i][j-dhalf] = m2[i][j];
				if (odd)
					t2[i][dhalf-1] = 0;
			}
			// decrement temp2 by H
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t2[i-dhalf][j-dhalf] -= m2[i][j];	
			
			strassen(t1, t2, p, dhalf);
			// update res with P1
			for (int i = 0; i < dhalf; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] = p[i][j-dhalf];
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] = p[i-dhalf][j-dhalf];
			
			
			
			// P2 = (A+B)H ===========================
			
			// A is already in t1
			// Add B to t1
			for (int i = 0; i < dhalf; i++)
				for (int j = dhalf; j < dim; j++)
					t1[i][j-dhalf] += m1[i][j];
			// set temp2 to H 
			// note if odd, rightmost column already padded
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t2[i-dhalf][j-dhalf] = m2[i][j];
			if (odd)
			{
				for (int j = dhalf; j < dim; j++)
					t2[dhalf-1][j-dhalf] = 0;
			}
			
			strassen(t1, t2, p, dhalf);
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] = -p[i][j];
			for (int i = 0; i < dhalf; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] += p[i][j-dhalf];
			
			
			// P3 = (C+D)E ===========================
			// copy C into t1
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					t1[i-dhalf][j] = m1[i][j];
			if (odd)
			{
				for (int j = 0; j < dhalf; j++)
					t1[dhalf-1][j] = 0;
			}
			// add D to t1
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t1[i-dhalf][j-dhalf] += m1[i][j];
			// copy E to t2
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t2[i][j] = m2[i][j];
			
			strassen(t1, t2, p, dhalf);
			// update result matrix with p3
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] = p[i-dhalf][j];
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] -= p[i-dhalf][j-dhalf];
			
			
			
			
			// P4 = D(G-E) ===========================
			// restore t1 to D
			// note if odd, bottom row already padded
			for (int i = dhalf; i < dim; i++)
			{
				for (int j = dhalf; j < dim; j++)
					t1[i-dhalf][j-dhalf] = m1[i][j];
				if (odd)
					t1[i-dhalf][dhalf-1] = 0;
			}
			// set t2 to G-E. E is already in t2.
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					t2[i-dhalf][j] = m2[i][j] - t2[i-dhalf][j];
			if (odd)
			{
				for (int j = 0; j < dhalf; j++)
					t2[dhalf-1][j] *= -1;
			}
			
			strassen(t1, t2, p, dhalf);
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] += p[i][j];
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] += p[i-dhalf][j];
			
			
		
			// P5 = (A+D)(E+H) =======================
			// set t1 to (A+D) by adding A to t1 as t1 = D already
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t1[i][j] += m1[i][j];
			// set t2 to E
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t2[i][j] = m2[i][j];
			// add H to t2
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t2[i-dhalf][j-dhalf] += m2[i][j];
			
			strassen(t1, t2, p, dhalf);
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] += p[i][j];
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] += p[i-dhalf][j-dhalf];
			
			
			// P6 = (B-D)(G+H) =======================
			// copy B into t1
			for (int i = 0; i < dhalf; i++)
			{
				for (int j = dhalf; j < dim; j++)
					t1[i][j-dhalf] = m1[i][j];
				if (odd)
					t1[i][dhalf-1] = 0;
			}
			// subtract D from t1
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t1[i-dhalf][j-dhalf] -= m1[i][j];
			// copy G into t2
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					t2[i-dhalf][j] = m2[i][j];
			if (odd)
			{
				for (int j = 0; j < dhalf; j++)
					t2[dhalf-1][j] = 0;
			}
			// add H to t2
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					t2[i-dhalf][j-dhalf] += m2[i][j];
		
			strassen(t1, t2, p, dhalf);
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					res[i][j] += p[i][j];
			
			// P7 = (A-C)(E+F) =======================
			// copy A into t1
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t1[i][j] = m1[i][j];
			// subtract C from t1
			for (int i = dhalf; i < dim; i++)
				for (int j = 0; j < dhalf; j++)
					t1[i-dhalf][j] -= m1[i][j];
			// copy E into t2
			for (int i = 0; i < dhalf; i++)
				for (int j = 0; j < dhalf; j++)
					t2[i][j] = m2[i][j];
			// subtract F from t2
			for (int i = 0; i < dhalf; i++)
				for (int j = dhalf; j < dim; j++)
					t2[i][j-dhalf] += m2[i][j];
			
			strassen(t1, t2, p, dhalf);
			for (int i = dhalf; i < dim; i++)
				for (int j = dhalf; j < dim; j++)
					res[i][j] -= p[i-dhalf][j-dhalf];
			
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
		System.out.println("");
	}
	
	public static boolean equals(int[][] m1, int[][] m2, int dim)
	{
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				if (m1[i][j] != m2[i][j])
					return false;
			}
		}
		return true;
	}
	
	public static void mult(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++) 
				res[i][j] = m1[i][0] * m2[0][j];		

		for (int k = 1; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] += (m1[i][k] * m2[k][j]);		
				}
	}
	
	public static void readMatrices(String filename, int[][] a, int[][] b) 
	{
		try 
		{
			Scanner sc = new Scanner(new File(filename));
			
			for (int i = 0; i < DIMENSION; i++)
				for (int j = 0; j < DIMENSION; j++) 
					a[i][j] = sc.nextInt();
				
			for (int i = 0; i < DIMENSION; i++)
				for (int j = 0; j < DIMENSION; j++) 
					b[i][j] = sc.nextInt();
			
			sc.close();
		} 
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
			return;
		}
	}

	public static void printDiag(int[][] m)
	{
		for (int i = 0; i < m.length; i++)
			System.out.println(m[i][i]);
	}
	
}
