/* -------------------
* FPAUspliceone.java
* -------------------
* This is the executable part of the FPAU project.  
* It uses logistic splice with historical fits to choose any carrying capacities.
*
*/
//package fpamodel;

import java.io.*;
//import java.util.StringTokenizer;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class FPAUspliceone
{
	// time span of the simulation [0,tspan]
	private int tspan = 2000;  //origin of time will be the year 1600 (i.e. 600 in our scale).  tspan is number of years to simulate after the year 1009.
	// External parameters
	private double alpha = 0.41;//0.45;
	private double beta = 0.001;
	private double gamma = 0.4;
	private double delta = 1.0;
	private double epsilon = 0;
	private double zeta = 0.9;
	/* Internal parameters for logistic inputs of the form
	 * x = K_x / (1 + exp(-r_x(t-t_Ix))) + x_0
	 * where x = p (population), c (consumption), and y (yield)
	 */
	//private double Kp,rp,tIp;
	//private double p0 = 6839652394.0; //population in 2009
	//private double dp0 = 78600758.36;
	private double Kp = 10486500000.0;
	private double rp = 0.031902;
	private double tIp = 13.2986+980;
	private double p0 = 310000000;
	
	//private double Kc,rc,tIc;
	private double c0 = 1659.959446; //consumption in 2009
	private double dc0 = 8.924201503;
	//private double Kc = 1940.6;//1609.8;//1925.1;
	//private double rc = 0.018671;//0.023148;//0.019399;;
	//private double tIc = 10.8366+985;//-12.7237+985;//-17.5155+980;
	//private double c0 = 571;//530;//328.0;

	//private double Ky,ry,tIy;
	private double y0 = 2274.133385; //yield value in 2009
	private double dy0 = 30.89071766;
	//private double Ky = 3390.6;//9162.8;//5050.5;//3547.6;//3447.6;
	//private double ry = 0.038934;//0.03637;//0.036148;//0.03719;//0.037202;
	//private double tIy = 10.7177+985;//5.8748+985;//12.4573+985;//11.6055+985;//7.3933+985;
	//private double y0 = 150;//85;//80.0;
	
	// the carrying capacities you get to choose.
	private double Ky = 3390.6+1500.0;
	private double Kc = 1940.6+571.0;
	
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
	private double[] c = new double[tspan+410];//remember origin of time is 600, so 410 corresponds to end of data, i.e. 2009
	private double[] p = new double[tspan+410];
	private double[] y = new double[tspan+410];
	// Initial conditions
	private double AF0 = 0.04306431889929523*T;//alpha*c0*p0/y0; //0.5 is fraction of demand coming from forests
	private double AP0 = 0.06197060524532735*T;//(1-alpha)*c0*p0/y0; //0.5 is fraction of demand coming from pasture
	private double s = 0.061;
	private double U0 = 0.0016520872299080305;//s*p0; //0.061 is ha/person using 3% urban area at 6.43 billion people in 2004
	private double BF0 = 0.0;
	private double BP0 = 0.0;
	private double F0 = 0.5254488025937881*T;//0.57*T - AF0 - zeta*U0;//0.6183*T - AF0;//0.637*T - AF0; //0.637 comes from biome analysis of Olson's map.
	private double P0 = 0.36786418603168225*T;//0.43*T - AP0 - (1-zeta)*U0;//0.3817*T - AF0;//0.363*T - AP0;
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
		FPAUspliceone fpa = new FPAUspliceone();
		fpa.initializeSimulation();
		fpa.runSimulation();
		fpa.outputResults();
		//fpa.outputResults();
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
			
			if (F[t+1]<0)
			{
				F[t+1]=0;
				AF[t+1]=AF[t+1]-alpha*Rplus*Theta(Rplus);
			}
			if (P[t+1]<0)
			{
				P[t+1]=0;
				AP[t+1]=AP[t+1]-alpha*Rplus*Theta(Rplus);
			}
			if (AF[t+1]<0)
			{
				AF[t+1]=0;
			}
			if (AP[t+1]<0)
			{
				AP[t+1]=0;
			}
		}
	}
	
	private void initializeSimulation()
	{
		//In principle, we could save a factor of n performance by updating these inputs at simulation time.  
		//But for now, it is conceptually simpler for me to do this initialization separately, although it is grossly inefficient.
		
		// load up the initial conditions
		F[0] = F0;
		P[0] = P0;
		AF[0] = AF0;
		AP[0] = AP0;
		BF[0] = BF0;
		BP[0] = BP0;
		U[0] = U0;
		
		// load up the input functions
		// find splicing params for yield function
		
		double tshift = 1009;
		double ry = dy0*Ky/y0/(Ky-y0);
		double tIy = (dy0*Ky*tshift+y0*Ky*Math.log((Ky-y0)/y0)-y0*y0*Math.log((Ky-y0)/y0))/dy0/Ky;
		// find splicing params for consumption function
		double rc = dc0*Kc/c0/(Kc-c0);
		double tIc = (dc0*Kc*tshift+c0*Kc*Math.log((Kc-c0)/c0)-c0*c0*Math.log((Kc-c0)/c0))/dc0/Kc;
		// initialize the input functions c and y
		
		System.out.println("ry: " + ry + "  tIy: " + tIy);
		System.out.println("rc: " + rc + "  tIc: " + tIc);
		
		//historical fitting parameters
		double Kch = 1940.6;//1609.8;//1925.1;
		double rch = 0.018671;//0.023148;//0.019399;;
		double tIch = 10.8366+985;//-12.7237+985;//-17.5155+980;
		double c0h = 571;//530;//328.0;

		double Kyh = 3390.6;//9162.8;//5050.5;//3547.6;//3447.6;
		double ryh = 0.038934;//0.03637;//0.036148;//0.03719;//0.037202;
		double tIyh = 10.7177+985;//5.8748+985;//12.4573+985;//11.6055+985;//7.3933+985;
		double y0h = 150;//85;//80.0;
		
		//XYSeries cseries = new XYSeries("C");
		//XYSeries yseries = new XYSeries("Y");
		//XYSeries pseries = new XYSeries("P");
		
		//XYSeries CPseries = new XYSeries("CP");
		//XYSeries Yseries = new XYSeries("YT");
		for (int t=0; t<410+tspan; t++)//remember origin of time is 600, so 410 corresponds to end of data, i.e. 2009
		{
			if (t<410)
			{
				c[t] = Kch / (1 + Math.exp(-rch*(t+600-tIch))) + c0h;
				//p[t] = Kp / (1 + Math.exp(-rp*(t+600-tIp))) + p0;
				y[t] = Kyh / (1 + Math.exp(-ryh*(t+600-tIyh))) + y0h;
			}
			else
			{
				//c[t] = 8.922810900*(t-410) + 1668.852953; //linear
				//y[t] = 30.89071766*(t-410) + 2304.868870; //linear
				c[t] = Kc / (1 + Math.exp(-rc*(t+tshift-410+1-tIc)));
				////p[t] = Kp / (1 + Math.exp(-rp*(t+1+tshift-tIp)))+p0;//not spliced
				y[t] = Ky / (1 + Math.exp(-ry*(t+tshift-410+1-tIy)));
			}
			p[t] = Kp / (1 + Math.exp(-rp*(t+600-tIp))) + p0;
			
			//cseries.add(t,c[t]/Kc);
			//yseries.add(t,y[t]/Ky);
			//pseries.add(t,p[t]/Kp);
			
			//CPseries.add(t,c[t]*p[t]/T);
			//Yseries.add(t,y[t]);
				
		}
		//XYSeriesCollection cpydataset = new XYSeriesCollection();
		/////cpydataset.addSeries(cseries);
		//XYSeriesCollection ydataset = new XYSeriesCollection();
		//cpydataset.addSeries(yseries);
		//XYSeriesCollection pdataset = new XYSeriesCollection();
		//////cpydataset.addSeries(pseries);	
		
		//ScatterPlot plotcpy = new ScatterPlot("Input Functions",cpydataset,800,500);
		//plotcpy.pack();
		//RefineryUtilities.centerFrameOnScreen(plotcpy);
		//plotcpy.setVisible(true);
		
		//XYSeriesCollection CPYdataset = new XYSeriesCollection();
		//CPYdataset.addSeries(CPseries);
		//CPYdataset.addSeries(Yseries);
		
		//ScatterPlot plotCPY = new ScatterPlot("CPYT",CPYdataset,800,500);
		//plotCPY.pack();
		//RefineryUtilities.centerFrameOnScreen(plotCPY);
		//plotCPY.setVisible(true);
		
		
		/*try
		{
			// initialize the file writer and dataset
			dataWriter = new PrintWriter(
					new BufferedWriter(
							new FileWriter("FPAFAPDspliceone.dat")));
		}
		catch (IOException e1)
		{
			System.out.println("Error creating FPAFAPDspliceone.dat");
		}*/
	}
	
	private void outputResults()
	{
		try
		{
			// initialize the file writer and dataset
			String pwd = System.getProperty("user.dir");
			dataWriter = new PrintWriter(
					new BufferedWriter(
							new FileWriter(pwd+"/data/output.dat")));
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
			
			//XYSeries CPseries = new XYSeries("dCP");
			//XYSeries YAseries = new XYSeries("dYA");
			//XYSeries CPYseries = new XYSeries("dCPY");
			
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
				//System.out.println("Writing " + t);
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
				
				//CPseries.add(t,(c[t+1]*p[t+1]-c[t]*p[t])/1000000000);
				//YAseries.add(t,y[t+1]-y[t]);
				//CPYseries.add(t,c[t]*p[t]/y[t]);
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
			
			//XYSeriesCollection CPYdataset = new XYSeriesCollection();
			//CPYdataset.addSeries(CPYseries);
			//CPYdataset.addSeries(YAseries);
			
			//ScatterPlot plotCPY = new ScatterPlot("CPY",CPYdataset,800,500);
			//plotCPY.pack();
			//RefineryUtilities.centerFrameOnScreen(plotCPY);
			//plotCPY.setVisible(true);
		}
		catch (IOException e1)
		{
			System.out.println("Error creating output.dat");
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
