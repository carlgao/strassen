import java.util.Random;
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class strassen 
{
	public static final int THRESHOLD = 8;
	
	private static int[][] a, b, cnorm, cstr; // TODO what do the last two do? also, note to self: p_i order optimization

	private static int d;

	public static void main (String[] args)
	{
		System.out.println(args.length);
		if (args.length != 3) {
			System.out.println("Usage: ./strassen 0 dimension inputfile");
			return;
		}
		int dim = Integer.parseInt(args[1]);
		dim += dim & 1;
		
		a = randMatrix(dim, -2, 2);
		b = randMatrix(dim, -2, 2);
		int[][] res = new int[dim][dim];
		strassen(a, b, res, 0, 0, 0, 0, dim);
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				System.out.print(res[i][j] + " ");
			}
			System.out.println();
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
	
	public static void getArrays(String filename)
	{
		Scanner sc = null;
		try {
			sc = new Scanner(new File(filename));
		} catch (FileNotFoundException e) {
			e.printStackTrace();  
		}

		// populate a
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
				a[i][j] = sc.nextInt();

		// populate b
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++)
				b[i][j] = sc.nextInt();

		sc.close();
	}

	// cache-efficient standard multiplication algorithm
	public static void mult(int[][] m1, int[][] m2, int[][] res, int dim)
	{
		for (int k = 0; k < dim; k++)
			for (int i = 0; i < dim; i++)
				for (int j = 0; j < dim; j++)
					res[i][j] = res[i][j] + (m1[i][k] * m2[k][j]);		
	}


	// nub memory usage dynamic padding style yay
	// RIGHT NOW THIS DOES NOT WORK CRY
	public static void strassen(int[][] m1, int[][] m2, int[][] res, int or1, int oc1, int or2, int oc2, int dim)
	{
		if (dim < THRESHOLD)
			mult(m1, m2, res, dim);
		else
		{
			int dm = (dim+1)>>1; // account for padding
			int diff = dim - dm;


			int[][] temp1 = new int[dm][dm];
			int[][] temp2 = new int[dm][dm];


			// P1
			int[][] p1 = new int[dm][dm];		
			// copy f-h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp2[i][j] = m2[i][j+dm] - m2[i+dm][j+dm];
				}
			}
			strassen(m1, temp2, p1, or1, oc1, 0, 0, dm);


			// P2
			int[][] p2 = new int[dm][dm];		
			// copy a+b into temp1
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i][j] + m1[i][j+dm];
				}
			}
			strassen(temp1, m2, p2, 0, 0, or2 + dm, oc2 + dm, dm);


			// P3
			int[][] p3 = new int[dm][dm];		
			// copy c+d into temp1
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i+dm][j] + m1[i+dm][j+dm];
				}
			}
			strassen(temp1, m2, p3, 0, 0, or2, oc2, dm);

			// P4
			int[][] p4 = new int[dm][dm];		
			// copy g-e into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp2[i][j] = m2[i+dm][j] - m2[i][j];
				}
			}
			strassen(m1, temp2, p4, or1+dm, oc1+dm, 0, 0, dm);

			// P5
			int[][] p5 = new int[dm][dm];
			// copy a+d into temp1, e+h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i][j] + m1[i+dm][j+dm];
					temp2[i][j] = m2[i][j] + m2[i+dm][j+dm];
				}
			}		
			strassen(temp1, temp2, p5, 0, 0, 0, 0, dm);

			// P6
			int[][] p6 = new int[dm][dm];
			// copy b-d into temp1, g+h into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i][j+dm] - m1[i+dm][j+dm];
					temp2[i][j] = m2[i+dm][j] + m2[i+dm][j+dm];
				}
			}		
			strassen(temp1, temp2, p6, 0, 0, 0, 0, dm);

			// P7
			int[][] p7 = new int[dm][dm];
			// copy a-c into temp1, e+f into temp2
			for (int i = 0; i < diff; i++)
			{
				for (int j = 0; j < diff; j++)
				{
					temp1[i][j] = m1[i][j] - m1[i+dm][j];
					temp2[i][j] = m2[i][j] + m2[i][j+dm];
				}
			}		
			strassen(temp1, temp2, p7, 0, 0, 0, 0, dm);

			// calculate result matrix
			for (int i = 0; i < dm; i++)
			{
				for (int j = 0; j < dm; j++)
				{
					res[i][j] = p5[i][j] + p4[i][j] - p2[i][j] + p6[i][j];
					res[i][j+dm] = p1[i][j] + p2[i][j];
					res[i+dm][j] = p3[i][j] + p4[i][j];
					res[i+dm][j+dm] = p5[i][j] + p1[i][j] - p3[i][j] - p7[i][j];
				}
			}
		}
	}
}
