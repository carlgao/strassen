import java.util.Random;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;

public class StrassenTimed
{
	public static final int THRESHOLD = 2; // must be >= 2
	public static final int TRIALS = 100;
	public static Map<Integer, int[][]> ps = new HashMap<Integer, int[][]>();
	public static Map<Integer, int[][]> temp1s = new HashMap<Integer, int[][]>();
	public static Map<Integer, int[][]> temp2s = new HashMap<Integer, int[][]>();
	
	/*
	 * Initializes temporary storage
	 */
	public void initialize()
	{
		
	}
	public static void main (String[] args)
	{
		boolean STRASSEN = true;
		int DIMENSION = 2; //Integer.parseInt(args[0]);
		
		int d2 = DIMENSION + (DIMENSION & 1);
		
		int dim = d2;
		
		// initialize storage for all dimensions of strassen's recursion.
		while (dim > THRESHOLD) 
		{
			dim = dim + (dim & 1);
			ps.put(dim, new int[dim][dim]);
			temp1s.put(dim, new int[dim][dim]);
			temp2s.put(dim, new int[dim][dim]);
			dim >>= 1;
		}
		// add base case of recursion
		ps.put(dim, new int[dim][dim]);
		temp1s.put(dim, new int[dim][dim]);
		temp2s.put(dim, new int[dim][dim]);
		
		// check keys 
		Iterator<Integer> iterator = ps.keySet().iterator();
		System.out.print("keys: ");
		while (iterator.hasNext())
			System.out.print((Integer)iterator.next()+" ");
		System.out.println("");

		int[][] a = randMatrix(d2, 1, 1);
		int[][] b = randMatrix(d2, 0, 1);

		// if odd, pad the matrices with 0s
		if (1 == (DIMENSION & 1)) 
		{
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
		}

		printMatrix(a);
		printMatrix(b);

		// TIMING SCRIPT
		/*
		
		long t = 0;
		long ttemp = -1;
		
		// discard first 10 results because the program seems to be slow when "warming up" initially
		for (int i = 0; i < 10; i++)
		{
			int[][] c = new int[d2][d2];
			int[][] d = new int[d2][d2];
			
			if (STRASSEN)
			{
				ttemp = System.nanoTime();
				strassen(a, b, d, DIMENSION);
				ttemp = System.nanoTime() - ttemp;
				//printMatrix(d);
			}
			else 
			{
				ttemp = System.nanoTime();
				mult(a, b, c, DIMENSION);
				ttemp = System.nanoTime() - ttemp;
				//printMatrix(c);
				//System.out.println("");
			}
		}
		
		for (int i = 0; i < TRIALS; i++)
		{
			int[][] c = new int[d2][d2];
			int[][] d = new int[d2][d2];
			
			if (STRASSEN)
			{
				ttemp = System.nanoTime();
				strassen(a, b, d, DIMENSION);
				ttemp = System.nanoTime() - ttemp;
				//printMatrix(d);
			}
			else 
			{
				ttemp = System.nanoTime();
				mult(a, b, c, DIMENSION);
				ttemp = System.nanoTime() - ttemp;
				//printMatrix(c);
				//System.out.println("");
			}
			t += ttemp;
		}
		t /= TRIALS;

		if (STRASSEN)
		{
			System.out.println("STRASSEN");
			System.out.println("Threshold: "+THRESHOLD);
		}
		else
		{
			System.out.println("TRADITIONAL");
		}
		
		System.out.println("Dimension: "+DIMENSION);
		System.out.println("Average computation time: "+ t + " ns");
		*/
		
		int[][] c = new int[d2][d2];
		int[][] d = new int[d2][d2];
		
		mult(a, b, c, DIMENSION);
		strassen(a, b, d, DIMENSION);
		
		System.out.println("c");
		printMatrix(c);
		System.out.println("d");
		printMatrix(d);

		for (int i = 0; i < DIMENSION; i++)
		{
			for (int j = 0; j < DIMENSION; j++)
			{
				if (c[i][j] != d[i][j])
				{
					System.out.println("Different at "+i+" "+j);
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
	 */
	
	public static void strassen(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		if (dim < THRESHOLD)
			strasMult(m1, m2, res, dim);
		
		else
		{
			// account for padding
			int dm = (dim+1) >> 1; 
			int dm2 = dm + (dm & 1);

			// temporary storage for factors of P matrices
			int[][] temp1 = temp1s.get(dm2);
			int[][] temp2 = temp2s.get(dm2);
			
			// temporary storage for P matrices
			int[][] p = ps.get(dm2);	
						
			// need to reset padding in case temps were used already
			if ((dm & 1) == 1)
			{
				for (int i = 0; i < dm; i++)
				{
					temp1[i][dm] = 0;
					temp1[dm][i] = 0;
					temp2[i][dm] = 0;
					temp2[dm][i] = 0;
					p[i][dm] = 0;
					p[dm][i] = 0;
				}
			}

			// P1 = A(F-H)
			// copy A into temp1 & F-H into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					//System.out.println("P1: "+"m2["+(i)+","+(j+dm)+"]="+m2[i][j+dm]+"  m2["+(i+dm)+","+(j+dm)+"]="+m2[i+dm][j+dm]);
					temp1[i][j] = m1[i][j];
					temp2[i][j] = m2[i][j+dm] - m2[i+dm][j+dm];
				}
			}
			
			strassen(temp1, temp2, p, dm);
			
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j+dm] = p[i][j];
					res[i+dm][j+dm] = p[i][j];
				}
			}
			System.out.println("P1");
			printMatrix(res);


			// P2 = (A+B)H
			// Copy A+B into temp1 & H into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j] + m1[i][j+dm];
					temp2[i][j] = m2[i+dm][j+dm];
				}
			}
			strassen(temp1, temp2, p, dm);
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j] = -p[i][j];
					res[i][j+dm] += p[i][j];
				}
			}
			System.out.println("P2");
			printMatrix(res);

			
			// P3 = (C+D)E
			// Copy C+D into temp1 & E into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i+dm][j] + m1[i+dm][j+dm];
					temp2[i][j] = m2[i][j];
				}
			}
			strassen(temp1, temp2, p, dm);
			
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i+dm][j] = p[i][j];
					res[i+dm][j+dm] -= p[i][j];
				}
			}
			
			// P4 = D(G-E)
			// Copy D into temp1 & G-E into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i+dm][j+dm];
					temp2[i][j] = m2[i+dm][j] - m2[i][j];
				}
			}
			strassen(temp1, temp2, p, dm);

			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j] -= p[i][j];				
					res[i+dm][j] += p[i][j];
				}
			}

			// P5 = (A+D)(E+H)
			// copy A+D into temp1 & E+H into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j] + m1[i+dm][j+dm];
					temp2[i][j] = m2[i][j] + m2[i+dm][j+dm];
				}
			}			
			strassen(temp1, temp2, p, dm);
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j] += p[i][j];	
					res[i+dm][j+dm] += p[i][j];
				}
			}
	
			// P6 = (B-D)(G+H)
			// Copy B-D into temp1 & G+H into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j+dm] - m1[i+dm][j+dm];
					temp2[i][j] = m2[i+dm][j] + m2[i+dm][j+dm];
				}
			}		
			strassen(temp1, temp2, p, dm);
	
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j] += p[i][j];
				}
			}
	
			// P7 = (A-C)(E+F)
			// Copy A-C into temp1 & E+F into temp2
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					temp1[i][j] = m1[i][j] - m1[i+dm][j];
					temp2[i][j] = m2[i][j] + m2[i][j+dm];
				}
			}		
			strassen(temp1, temp2, p, dm);
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i+dm][j+dm] -= p[i][j];
				}
			}
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

	/*
	 * Cache-efficient standard multiplication algorithm. Assumes result
	 * matrix res is initialized to all 0's.
	 */
	public static void mult(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		for (int k = 0; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] = res[i][j] + (m1[i][k] * m2[k][j]);		
				}
	}

	/*
	 * Standard multiplication algorithm, used by Strassen when n < 0.
	 * Does not require res to be all 0.
	 */
	public static void strasMult(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++) 
				res[i][j] = (m1[i][0] * m2[0][j]);		

		for (int k = 1; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++) {
					res[i][j] += (m1[i][k] * m2[k][j]);		
				}
	}
	
	/*
	 * Prints the diagonal of a square matrix m, one entry per line
	 */
	public static void printDiag(int[][] m)
	{
		for (int i = 0; i < m.length; i++)
			System.out.println(m[i][i]);
	}
}

