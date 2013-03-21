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
		if (args.length != 3)
		{
			System.out.println("Usage: ./strassen 0 dimension inputfile");
			return;
		}
		d = Integer.parseInt(args[1]);
		int d2 = d + (d & 1);
		
		int[][] a = new int[d2][d2];
		int[][] b = new int[d2][d2];
		getArrays(args[2], a, b);
		int[][] res = new int[d2][d2];
		strassen(a, b, res, 0, 0, 0, 0, d2);
		for (int i = 0; i < d; i++)
		{
			for (int j = 0; j < d; j++)
			{
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
	
	public static void getArrays(String filename, int[][] a, int[][] b) // why can't I access a and b without passing them in??
	{
		Scanner sc = null;
		try {
			sc = new Scanner(new File(filename));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}

		// populate a
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++) {
				a[i][j] = sc.nextInt();
			}

		// populate b
		for (int i = 0; i < d; i++)
			for (int j = 0; j < d; j++) {
				b[i][j] = sc.nextInt();
			}

		sc.close();
	}

	// cache-efficient standard multiplication algorithm
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


	// nub memory usage dynamic padding style yay
	// RIGHT NOW THIS DOES NOT WORK CRY
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
}

