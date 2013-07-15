/* -------------------
* FPAU.java
* -------------------
* This is the executable part of the FPAU project.  It loops through a range of carrying capacities and splices logistic curves to historical fit.
*
*/
package fpamodel;

import java.io.*;
//import java.util.StringTokenizer;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class FPAUphases
{
	// time span of the simulation [0,tspan]
	private int tspan = 2000;
	// External parameters
	private double alpha = 0.4;//0.45;
	private double beta = 0.01;
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
	private double AF0 = 0.1773529548330691*T;//alpha*c0*p0/y0; //0.5 is fraction of demand coming from forests
	private double AP0 = 0.2660294322496037*T;//(1-alpha)*c0*p0/y0; //0.5 is fraction of demand coming from pasture
	private double s = 0.061;
	private double U0 = 0.03705370716776464;//s*p0; //0.061 is ha/person using 3% urban area at 6.43 billion people in 2004
	private double BF0 = 0.0;
	private double BP0 = 0.0;
	private double F0 = 0.3592987087159442*T;//0.57*T - AF0 - zeta*U0;//0.6183*T - AF0;//0.637*T - AF0; //0.637 comes from biome analysis of Olson's map.
	private double P0 = 0.16026519703362008*T;//0.43*T - AP0 - (1-zeta)*U0;//0.3817*T - AF0;//0.363*T - AP0;
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
		FPAUphases fpa = new FPAUphases();
		fpa.initializeSimulation();
		fpa.runSimulation();
		//fpa.outputResults();
	}
	
	private void runSimulation()
	{	
		double Kymin = 2400;
		double Kymax = 10000;
		double Kcmin = 1700;
		double Kcmax = 6000;
		
		double ry,tIy,rc,tIc;
		
		double tshift = 1009;
		
		double iter = 0;
		
		for (double Ky=Kymin; Ky <= Kymax; Ky=Ky+50)
		{
			// find splicing params for yield function
			ry = dy0*Ky/y0/(Ky-y0);
			tIy = (dy0*Ky*tshift+y0*Ky*Math.log((Ky-y0)/y0)-y0*y0*Math.log((Ky-y0)/y0))/dy0/Ky;
			for (double Kc=Kcmin; Kc<=Kcmax; Kc=Kc+25)
			{
				iter++;//useless.
				System.out.println(Ky + "\t" + Kc);
				// find splicing params for consumption function
				rc = dc0*Kc/c0/(Kc-c0);
				tIc = (dc0*Kc*tshift+c0*Kc*Math.log((Kc-c0)/c0)-c0*c0*Math.log((Kc-c0)/c0))/dc0/Kc;
				// initialize the input functions c and y
				c[0] = Kc / (1 + Math.exp(-rc*(0+tshift-tIc)));
				p[0] = Kp / (1 + Math.exp(-rp*(0+tshift-tIp)));//not spliced
				y[0] = Ky / (1 + Math.exp(-ry*(0+tshift-tIy)));
				
				// load up the initial conditions
				F[0] = F0;
				P[0] = P0;
				AF[0] = AF0;
				AP[0] = AP0;
				BF[0] = BF0;
				BP[0] = BP0;
				U[0] = U0;
				for (int t=0; t<tspan-1; t++)
				{
						c[t+1] = Kc / (1 + Math.exp(-rc*(t+1+tshift-tIc)));
						p[t+1] = Kp / (1 + Math.exp(-rp*(t+1+tshift-tIp)))+p0;//not spliced
						y[t+1] = Ky / (1 + Math.exp(-ry*(t+1+tshift-tIy)));
					
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
					}
					if (P[t+1]<0)
					{
						P[t+1]=0;
					}
					if (AF[t+1]<0)
					{
						AF[t+1]=0;
					}
					if (AP[t+1]<0)
					{
						AP[t+1]=0;
					}
					//normalize previous result
					F[t] = F[t]/T;
					P[t] = P[t]/T;
					AF[t] = AF[t]/T;
					AP[t] = AP[t]/T;
					BF[t] = BF[t]/T;
					BP[t] = BP[t]/T;
					D[t] = (c[t]*p[t]/y[t]-(AF[t]+AP[t]))/T;
					U[t] = U[t]/T;
				}
				F[tspan-1] = F[tspan-1]/T;
				P[tspan-1] = P[tspan-1]/T;
				AF[tspan-1] = AF[tspan-1]/T;
				AP[tspan-1] = AP[tspan-1]/T;
				BF[tspan-1] = BF[tspan-1]/T;
				BP[tspan-1] = BP[tspan-1]/T;
				D[tspan-1] = (c[tspan-1]*p[tspan-1]/y[tspan-1]-(AF[tspan-1]+AP[tspan-1]))/T;
				U[tspan-1] = U[tspan-1]/T;
				getMetrics(Ky,Kc);
				if (iter==3070)
				{
					outputResults();
				}
			}
		}
		dataWriter.close();
	}
	
	private void getMetrics(double Ky, double Kc)
	{
		// This method could be integrated with runSimulation() to get the metrics at simulation time.  Done separately here for conceptual simplicity
		int numFTs=0;
		int numPTs=0;
		int numATs=0; //number of transitions
		int tFT=-1;
		int tPT=-1;
		int tAT=-1; //time of transitions
		int tFequib=-1; 
		int tPequib=-1;
		int tAequib=-1;  //times of onset of equilibrium.
		double FFT=2.0; 
		double PPT=2.0;
		double AAT=2.0; 
		double Fequib=2.0;
		double Pequib=2.0; 
		double Aequib=2.0; //land cover values at transitions and at equilibrium.
		int FequibTicker=0;
		int PequibTicker=0;
		int AequibTicker=0; //counts number of times slope is zero.
		
		double thisFslope,thisPslope,thisAslope;
		double lastFslope = F[1]-F[0];
		double lastPslope = P[1]-P[0];
		double lastAslope = (AF[1]+AP[1])-(AF[0]+AP[0]);
		
		int numfFTs=0;
		int numfPTs=0;
		int numfATs=0; //number of false transitions
		int tfFT=-1;
		int tfPT=-1;
		int tfAT=-1; //time of transitions
		double FfFT=2.0; 
		double PfPT=2.0;
		double AfAT=2.0;
		
		int phaseflag = -1; // 1 for no FT. 2 for FT. 3 for false FT. 4 for quasi-false FT.
		
		for (int t=1; t<tspan-1; t++)
		{
			thisFslope = F[t+1]-F[t];
			thisPslope = P[t+1]-P[t];
			thisAslope = (AF[t+1]+AP[t+1])-(AF[t]+AP[t]);
			//probe for transitions
			if (thisFslope > 0 && lastFslope <= 0)
			{
				numFTs++;
				if (numFTs==1)
				{
					tFT = t;
					FFT = F[t];
				}
				if (FequibTicker>0)//removes false equilibrium flags.
				{
					FequibTicker=0;
					tFequib=-1;
					Fequib=2.0;
				}
			}
			if (thisPslope > 0 && lastPslope <= 0)
			{
				numPTs++;
				if (numPTs==1)
				{
					tPT = t;
					PPT = P[t];
				}
				if (PequibTicker>0)//removes false equilibrium flags.
				{
					PequibTicker=0;
					tPequib=-1;
					Pequib=2.0;
				}
			}
			if (thisAslope < 0 && lastAslope >= 0)
			{
				numATs++;
				if (numATs==1)
				{
					tAT = t;
					AAT = AF[t]+AP[t];
				}
				if (AequibTicker>0)//removes false equilibrium flags.
				{
					AequibTicker=0;
					tAequib=-1;
					Aequib=2.0;
				}
			}
			// probe for false transitions
			if (thisFslope < 0 && lastFslope >= 0)
			{
				numfFTs++;
				if (numfFTs==1)
				{
					tfFT = t;
					FfFT = F[t];
				}
				if (FequibTicker>0)//removes false equilibrium flags.
				{
					FequibTicker=0;
					tFequib=-1;
					Fequib=2.0;
				}
			}
			if (thisPslope < 0 && lastPslope >= 0)
			{
				numfPTs++;
				if (numfPTs==1)
				{
					tfPT = t;
					PfPT = P[t];
				}
				if (PequibTicker>0)//removes false equilibrium flags.
				{
					PequibTicker=0;
					tPequib=-1;
					Pequib=2.0;
				}
			}
			if (thisAslope > 0 && lastAslope <= 0)
			{
				numfATs++;
				if (numfATs==1)
				{
					tfAT = t;
					AfAT = AF[t]+AP[t];
				}
				if (AequibTicker>0)//removes false equilibrium flags.
				{
					AequibTicker=0;
					tAequib=-1;
					Aequib=2.0;
				}
			}
			//probe for equilibria
			if (Math.abs(thisFslope)<Math.pow(10, -4))
			{
				FequibTicker++;
			}
			if (FequibTicker>100 && tFequib==-1)
			{
				tFequib = t-100;
				Fequib = F[t];
			}
			if (Math.abs(thisPslope)<Math.pow(10, -4))
			{
				PequibTicker++;
			}
			if (PequibTicker>100 && tPequib==-1)
			{
				tPequib = t-100;
				Pequib = P[t];
			}
			if (Math.abs(thisAslope)<Math.pow(10, -4))
			{
				AequibTicker++;
			}
			if (AequibTicker>100 && tAequib==-1)
			{
				tAequib = t-100;
				Aequib = AF[t]+AP[t];
			}
			if (tFequib!=-1 && tPequib!=-1 && tAequib!=-1)
			{
				t=tspan;
			}
			lastFslope=thisFslope;
			lastPslope=thisPslope;
			lastAslope=thisAslope;
		}
		// set the flags
		if (numFTs == 0)
		{
			phaseflag = 1;
			if (Fequib<0.00001)
			{
				phaseflag=0;
			}
		}
		if (numFTs>0 && numfFTs==0)
		{
			phaseflag = 2;
		}
		if (numFTs>0 && numfFTs>0 && Fequib-FFT<0.01)
		{
			phaseflag = 3;
		}
		if (numFTs>0 && numfFTs>0 && Fequib-FFT>=0.01)
		{
			phaseflag = 4;
		}
		// remove false false transitions that do not revert by more than 1%
		/*if (numfFTs>0 && Math.abs(FfFT-Fequib)<0.01)
		{
			numfFTs=0;
			FfFT=2.0;
			tfFT=-1;
		}
		if (numfPTs>0 && Math.abs(PfPT-Pequib)<0.01)
		{
			numfPTs=0;
			PfPT=2.0;
			tfPT=-1;
		}
		if (numfATs>0 && Math.abs(AfAT-Aequib)>0.01)
		{
			numfATs=0;
			AfAT=2.0;
			tfAT=-1;
		}
		// remove transitions that do not recover by at least 1%
		if (numFTs>0 && Math.abs(Fequib-FFT)<0.01)
		{
			numFTs=0;
			FFT=2.0;
			tFT=-1;
		}
		if (numPTs>0 && Math.abs(Pequib-PPT)<0.01)
		{
			numPTs=0;
			PPT=2.0;
			tPT=-1;
		}
		if (numATs>0 && Math.abs(Aequib-AAT)>0.01)
		{
			numATs=0;
			AAT=2.0;
			tAT=-1;
		}*/
		
		dataWriter.println(Ky+"\t"+Kc+"\t"+numFTs+"\t"+numPTs+"\t"+numATs+"\t"+tFT+"\t"+tPT+"\t"+tAT+"\t"+FFT+"\t"+PPT+"\t"+AAT+"\t"+tFequib+"\t"+tPequib+"\t"+tAequib+"\t"+Fequib+"\t"+Pequib+"\t"+Aequib+"\t"+numfFTs+"\t"+numfPTs+"\t"+numfATs+"\t"+tfFT+"\t"+tfPT+"\t"+tfAT+"\t"+FfFT+"\t"+PfPT+"\t"+AfAT+"\t"+phaseflag);
	
	}
	
	private void initializeSimulation()
	{
		try
		{
			// initialize the file writer and dataset
			dataWriter = new PrintWriter(
					new BufferedWriter(
							new FileWriter("FPAFAPDsplice.dat")));
		}
		catch (IOException e1)
		{
			System.out.println("Error creating FPAFAPDsplice.dat");
		}
	}
	
	private void outputResults()
	{
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
