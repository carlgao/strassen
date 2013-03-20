
public class Strassen 
{
	private int[][] a, b;
	private int d;
	
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
	public static int[][] mult()
	{
		int[][] c = new int[d][d];
		for (int k = 0; k < d; k++)
			for (int i = 0; i < d; i++)
				for (j = 0; j < d; j++)
					c[i][j] = c[i][j] + a[i][k]*b[k][j];		
		return c;
	}
}
