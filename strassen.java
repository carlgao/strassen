
public class Strassen 
{
	private int[][] a, b;
	private int n;
	
	public static void main (String[] args)
	{
		
	}
	
	public void getArrays(String filename)
	{
		Scanner sc = null;
	    try {
	        sc = new Scanner(new File(filename));
	    } catch (FileNotFoundException e) {
	        e.printStackTrace();  
	    }
	    
	    // populate a
	    for (int i = 0; i < n; i++)
	    	for (int j = 0; j < n; j++)
	    		a[i][j] = sc.nextInt();
	    
	    // populate b
	    for (int i = 0; i < n; i++)
	    	for (int j = 0; j < n; j++)
	    		b[i][j] = sc.nextInt();
	    
	    sc.close();
	}
	
	// cache-efficient standard multiplication algorithm
	public static int[][] mult()
	{
		int[][] c = new int[n][n];
		for (int k = 0; k < n; k++)
			for (int i = 0; i < n; i++)
				for (j = 0; j < n; j++)
					c[i][j] = c[i][j] + a[i][k]*b[k][j];		
		return c;
	}
}
