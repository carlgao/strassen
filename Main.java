import java.util.Random;


public class Main
{
	public static final int THRESHOLD = 55;
    
	public static void main (String[] args)
	{
		int DIMENSION = 140;//Integer.parseInt(args[0]);
		int d2 = DIMENSION + (DIMENSION & 1);
        
		int[][] a = randMatrix(d2, 0, 1);//{{1,0,0},{1,1,1},{0,0,1}};//
		int[][] b = randMatrix(d2, 1, 1);//{{0,0,1},{1,0,0},{1,1,1}};//randMatrix(2, 0, 1);
        
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
        
		long t1 = System.nanoTime();
		mult(a, b, c, DIMENSION);
		long t2 = System.nanoTime();
		strassen(a, b, d, DIMENSION);
		long t3 = System.nanoTime();
        
		if (t2-t1 > t3-t2)
			System.out.println("Strassen's is faster.");
		else
			System.out.println("Normal multiplication is faster.");
		System.out.println("Normal: "+(t2-t1)+"\t Strassen's: "+(t3-t2));
		
		/*
         System.out.println("A:");
         printMatrix(a);
         
         System.out.println("\nB:");
         printMatrix(b);
         
         System.out.println("\nStrassen's: ");
         printMatrix(d);
         
         System.out.println("\nActual: ");
         printMatrix(c);*/
        
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
			int[][] temp1 = new int[dm2][dm2];
			int[][] temp2 = new int[dm2][dm2];
            
			// temporary storage for P matrices
			int[][] p = new int[dm2][dm2];
            
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
            
            
			// P2 = (A+B)H
			// copy a+b into temp1 and h into temp2
			
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

