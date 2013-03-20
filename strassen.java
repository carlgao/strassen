
public class Strassen 
{
	private int[][] a, b;
	private int n;
	
	public static void main (String[] args)
	{
		
	}
	
	// cache-efficient standard multiplication algorithm
	public static int[][] mult()
	{
		int[][] c = new int[n][n];
		for (int k = 0; k < n; k++)
		{
			for (int i = 0; i < n; i++)
			{
				for (j = 0; j < n; j++)
				{
					c[i][j] = c[i][j] + a[i][k]*b[k][j];
				}
			}
		}		
		return c;
	}
}
