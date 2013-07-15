/* -------------------
* FPAU.java
* -------------------
* This is the executable part of the FPAU project.  It uses logistic fits to historical data.
*
*/
package fpamodel;

import java.io.*;
import java.util.StringTokenizer;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class FPAU 
{
	// time span of the simulation [0,tspan]
	private int tspan = 1500;
	// External parameters
	private double alpha = 0.41;//0.45;
	private double beta = 0.01;
	private double gamma = 0.535;
	private double delta = 1.0;
	private double epsilon = 0;
	private double zeta = 0.9;
	/* Internal parameters for logistic inputs of the form
	 * x = K_x / (1 + exp(-r_x(t-t_Ix))) + x_0
	 * where x = p (population), c (consumption), and y (yield)
	 */
	private double Kp = 10486500000.0;
	private double rp = 0.031902;
	private double tIp = 13.2986+980;
	private double p0 = 310000000;
	
	private double Kc = 1940.6;//1609.8;//1925.1;
	private double rc = 0.018671;//0.023148;//0.019399;;
	private double tIc = 10.8366+985;//-12.7237+985;//-17.5155+980;
	private double c0 = 571;//530;//328.0;

	private double Ky = 3390.6;//9162.8;//5050.5;//3547.6;//3447.6;
	private double ry = 0.038934;//0.03637;//0.036148;//0.03719;//0.037202;
	private double tIy = 10.7177+985;//5.8748+985;//12.4573+985;//11.6055+985;//7.3933+985;
	private double y0 = 150;//85;//80.0;
	
	private double T = 11259839521.0941;//9365600169.31;//9341604305.9;//8996886620.5;//13049972500.0; //total land area in ha excluding Antarctica
	
	/*
	 * Data arrays for state variables F (forest area), P (pasture), 
	 * AF (agricultural area derived from forest), AP (agricultural area derived from pasture),
	 * BF (abandoned land originally forested), BP (abandoned land originally pasture)
	 * D (demand, consumption*population/yield)
	 */
	private double[] F = new double[tspan];
	private double[] P = new double[tspan];
	private double[] AF = new double[tspan];
	private double[] AP = new double[tspan];
	private double[] BF = new double[tspan];
	private double[] BP  = new double[tspan];
	private double[] U  = new double[tspan];
	private double[] D = new double[tspan];
	/*
	 * Data arrays for input variables c (per capita annual consumption), 
	 * p (population), y (yield)
	 * Note that c*p/y = demand
	 */
	private double[] c = new double[tspan];
	private double[] p = new double[tspan];
	private double[] y = new double[tspan];
	// Initial conditions
	private double AF0 = alpha*c0*p0/y0; //0.5 is fraction of demand coming from forests
	private double AP0 = (1-alpha)*c0*p0/y0; //0.5 is fraction of demand coming from pasture
	private double s = 0.06;
	private double U0 = s*p0; //0.061 is ha/person using 3% urban area at 6.43 billion people in 2004
	private double BF0 = 0.0;
	private double BP0 = 0.0;
	private double F0 = 0.57*T - AF0 - zeta*U0;//0.6183*T - AF0;//0.637*T - AF0; //0.637 comes from biome analysis of Olson's map.
	private double P0 = 0.43*T - AP0 - (1-zeta)*U0;//0.3817*T - AF0;//0.363*T - AP0;
	// Object used for writing data to file
	private PrintWriter dataWriter;
	/*
	 * Objects used as input to ScatterPlot for generating graphs.
	 * Suffixes correspond to float arrays F,P,AF,AP,D
	 */
	private XYSeries seriesF,seriesP,seriesAF,seriesAP,seriesA,seriesBF,seriesBP,seriesB,seriesU,seriesD;
	private XYSeriesCollection dataset;
	
	public static void main(String[] args) 
	{
		FPAU fpa = new FPAU();
		fpa.initializeSimulation();
		fpa.runSimulation();
		fpa.outputResults();
	}
	
	private void runSimulation()
	{
		for (int t=0; t<tspan-1; t++)
		{
			double Rplus = c[t+1]*p[t+1]/y[t+1]-AF[t]-AP[t];
			double Rminus = -Rplus;
			
			F[t+1] = F[t] - alpha*Rplus*Theta(Rplus) + beta*BF[t] - zeta*s*(p[t+1]-p[t]);
			P[t+1] = P[t] - (1-alpha)*Rplus*Theta(Rplus) + delta*BP[t] - (1-zeta)*s*(p[t+1]-p[t]);
			AF[t+1] = AF[t] + alpha*Rplus*Theta(Rplus) - gamma*Rminus*Theta(Rminus) - epsilon*AF[t];
			AP[t+1] = AP[t] + (1-alpha)*Rplus*Theta(Rplus) - (1-gamma)*Rminus*Theta(Rminus) - epsilon*AP[t];
			BF[t+1] = BF[t] + gamma*Rminus*Theta(Rminus) - beta*BF[t] + epsilon*AF[t];
			BP[t+1] = BP[t] + (1-gamma)*Rminus*Theta(Rminus) - delta*BP[t] + epsilon*AP[t];
			U[t+1] = U[t] + s*(p[t+1]-p[t]);
		}
	}
	
	private void initializeSimulation()
	{
		// populate the input function p with data where data exists
		int line = 0;
		try
		{
			FileInputStream fstream = new FileInputStream("population_hermite.dat");
			DataInputStream in = new DataInputStream(fstream);
			BufferedReader br = new BufferedReader(new InputStreamReader(in));
			String strLine;
			StringTokenizer st;
			//Read File Line By Line
			while ((strLine = br.readLine()) != null && line<tspan)   
			{
				st = new StringTokenizer(strLine, "\t");
				st.nextToken();
				double val = (Kp+p0)*Double.parseDouble(st.nextToken());
				p[line] = val;
				line++;
			}
			//Close the input stream
			in.close();
		}
		catch (IOException ioe)
		{
			System.out.println("Error reading population data file.");
		}
		// populate the input functions c and y
		for (int t=0; t<tspan; t++)
		{
			c[t] = Kc / (1 + Math.exp(-rc*(t-tIc))) + c0;
			if (true)//(t>line-1) //uncomment to use historical population data as input
			{
				p[t] = Kp / (1 + Math.exp(-rp*(t-tIp))) + p0;
			}
			y[t] = Ky / (1 + Math.exp(-ry*(t-tIy))) + y0;
			double diff = ((p[t])-(Kp / (1 + Math.exp(-rp*(t-tIp))) + p0))/(Kp / (1 + Math.exp(-rp*(t-tIp))) + p0);
			System.out.println(t + "  " + p[t] + "  " + (Kp / (1 + Math.exp(-rp*(t-tIp))) + p0) + "  " + diff);
		}
		// load up the initial conditions
		F[0] = F0;
		P[0] = P0;
		AF[0] = AF0;
		AP[0] = AP0;
		BF[0] = BF0;
		BP[0] = BP0;
		U[0] = U0;
	}
	
	private void outputResults()
	{
		try
		{
			// initialize the file writer and dataset
			dataWriter = new PrintWriter(
					new BufferedWriter(
							new FileWriter("FPAFAPD.dat")));
			seriesF = new XYSeries("F");
			seriesP = new XYSeries("P");
			seriesAF = new XYSeries("AF");
			seriesAP = new XYSeries("AP");
			seriesA = new XYSeries("A");
			seriesBF = new XYSeries("BF");
			seriesBP = new XYSeries("BP");
			seriesB = new XYSeries("B");
			seriesU = new XYSeries("U");
			seriesD = new XYSeries("D");
			// write data to file and load XYDataset for plot
			for (int t=0; t<tspan; t++)
			{
				// normalize results with respect to total area T
				F[t] = F[t]/T;
				P[t] = P[t]/T;
				AF[t] = AF[t]/T;
				AP[t] = AP[t]/T;
				BF[t] = BF[t]/T;
				BP[t] = BP[t]/T;
				D[t] = (c[t]*p[t]/y[t]-(AF[t]+AP[t]))/T;
				U[t] = U[t]/T;
				dataWriter.println(t+"\t"+F[t]+"\t"+P[t]+"\t"+(AF[t]+AP[t])+"\t"+AF[t]+"\t"+AP[t]+"\t"+(BF[t]+BP[t])+"\t"+BF[t]+"\t"+BP[t]+"\t"+U[t]+"\t"+D[t]);
				seriesF.add(t,F[t]);
				seriesP.add(t,P[t]);
				seriesA.add(t,(AF[t]+AP[t]));
				seriesB.add(t,(BF[t]+BP[t]));
				seriesAF.add(t,AF[t]);
				seriesAP.add(t,AP[t]);
				seriesBF.add(t,BF[t]);
				seriesBP.add(t,BP[t]);
				seriesU.add(t,U[t]);
				seriesD.add(t,D[t]);
			}
			dataWriter.close();
			
			dataset = new XYSeriesCollection();
			dataset.addSeries(seriesF);
			dataset.addSeries(seriesP);
			dataset.addSeries(seriesA);
			dataset.addSeries(seriesB);
			dataset.addSeries(seriesU);
			dataset.addSeries(seriesAF);
			dataset.addSeries(seriesAP);
			dataset.addSeries(seriesBF);
			dataset.addSeries(seriesBP);
			//dataset.addSeries(seriesD);
			
			// generate the plot
			//XYDataset dataset = createDataset();
			ScatterPlot plot = new ScatterPlot("Land areas",dataset,800,500);
			plot.pack();
			RefineryUtilities.centerFrameOnScreen(plot);
			plot.setVisible(true);
		}
		catch (IOException e1)
		{
			System.out.println("Error creating FPAFAPD.dat");
		}
	}
	
	private int Theta(double x)
	{
		if (x>0)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}
