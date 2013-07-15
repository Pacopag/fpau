/* -------------------
* FPAU.java
* -------------------
* This is the executable part of the FPAU project.  It uses logistic fits to historical data.
*
*/

import java.io.*;
import java.util.StringTokenizer;
import java.util.Hashtable;
import java.util.Enumeration;

import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RefineryUtilities;

public class FPAU 
{
	private Hashtable<String,Double> params;
	
	/*
	 * Data arrays for state variables F (forest area), P (pasture), 
	 * AF (agricultural area derived from forest), AP (agricultural area derived from pasture),
	 * BF (abandoned land originally forested), BP (abandoned land originally pasture)
	 * D (demand, consumption*population/yield)
	 */
	private double[] F;
	private double[] P;
	private double[] AF;
	private double[] AP;
	private double[] BF;
	private double[] BP;
	private double[] U;
	private double[] D;
	/*
	 * Data arrays for input variables c (per capita annual consumption), 
	 * p (population), y (yield)
	 * Note that c*p/y = demand
	 */
	private double[] c;
	private double[] p;
	private double[] y;
	// Initial conditions
	private double AF0; 
	private double AP0;
	private double U0;
	private double BF0;
	private double BP0;
	private double F0;
	private double P0;
	
	private double T; //Total land area
	private int tspan; //Duration of the simulation
	
	/*
	 * Objects used as input to ScatterPlot for generating graphs.
	 * Suffixes correspond to float arrays F,P,AF,AP,D
	 */
	private XYSeries seriesF,seriesP,seriesAF,seriesAP,seriesA,seriesBF,seriesBP,seriesB,seriesU,seriesD;
	private XYSeriesCollection dataset;
	
	public static void main (String[] args)
	{
		try
		{
			FPAU fpa = new FPAU();
			fpa.run(args);
		}
		catch (Exception e)
		{
			e.printStackTrace ();
		}
	}
	
	public void run(String[] args) throws Exception
	{	
		if (args.length==0) {	
			initializeParams();
			initializeInputFns(false,0.0,0.0);
			runSimulation();
			outputResults("FPAU.dat");
		}
		else if (args.length==3 && args[0].equals("splice")) {	
			initializeParams();
			initializeInputFns(true,Double.parseDouble(args[0]),
				Double.parseDouble(args[1]));
			runSimulation();
			outputResults("FPAU_"+args[0]+"_"+args[1]+".dat");
		}
		else if (args.length==1 && args[0].equals("phases"))
		{
			initializeParams();
			initializeInputFns(false,0.0,0.0);
			runPhaseAnalysis();
		}
		else if (args.length==3 && args[0].equals("sensitivity"))
		{
			initializeParams();
			initializeInputFns(false,0.0,0.0);
			runSensitivityAnalysis(Double.parseDouble(args[1]),
				Integer.parseInt(args[2]));
		}
	}
	
	private void initializeParams()
	{
		params = new Hashtable<String,Double>();
		params.put("alpha",0.41);
		params.put("beta",0.01);
		params.put("gamma",0.4);
		params.put("delta",1.0);
		params.put("epsilon",0.0);
		params.put("zeta",0.9);
		params.put("s",0.06);
		params.put("Kp",10486500000.0);
		params.put("rp",0.031902);
		params.put("tIp",13.2986+980);
		params.put("p0",310000000.0);
		params.put("Kc",1940.6);
		params.put("rc",0.018671);
		params.put("tIc",10.8366+985);
		params.put("c0",571.0);
		params.put("Ky",3390.6);
		params.put("ry",0.038934);
		params.put("tIy",10.7177+985);
		params.put("y0",150.0);
		
		tspan = 2500;//(int)(P("tspan"));
		T = 11259839521.0941;
		
		F = new double[tspan];
		P = new double[tspan];
		AF = new double[tspan];
		AP = new double[tspan];
		BF = new double[tspan];
		BP  = new double[tspan];
		U  = new double[tspan];
		D = new double[tspan];
		/*
		 * Data arrays for input variables c (per capita annual consumption), 
		 * p (population), y (yield)
		 * Note that c*p/y = demand
		 */
		c = new double[tspan];
		p = new double[tspan];
		y = new double[tspan];
		// Initial conditions
		AF0 = P("alpha")*P("c0")*P("p0")/P("y0"); //0.5 is fraction of demand coming from forests
		AP0 = (1-P("alpha"))*P("c0")*P("p0")/P("y0"); //0.5 is fraction of demand coming from pasture
		U0 = P("s")*P("p0"); //0.061 is ha/person using 3% urban area at 6.43 billion people in 2004
		BF0 = 0.0;
		BP0 = 0.0;
		F0 = 0.57*T - AF0 - P("zeta")*U0;//0.6183*T - AF0;//0.637*T - AF0; //0.637 comes from biome analysis of Olson's map.
		P0 = 0.43*T - AP0 - (1-P("zeta"))*U0;//0.3817*T - AF0;//0.363*T - AP0;
	}
	
	private void initializeInputFns(boolean splice, double Ky, double Kc)
	{
		// load up the initial conditions
		F[0] = F0;
		P[0] = P0;
		AF[0] = AF0;
		AP[0] = AP0;
		BF[0] = BF0;
		BP[0] = BP0;
		U[0] = U0;
		
		if (!splice)
		{
			// populate the input functions c and y
			for (int t=0; t<tspan; t++)
			{
				c[t] = P("Kc") / (1 + Math.exp(-P("rc")*(t-P("tIc")))) + P("c0");
				p[t] = P("Kp") / (1 + Math.exp(-P("rp")*(t-P("tIp")))) + P("p0");
				y[t] = P("Ky") / (1 + Math.exp(-P("ry")*(t-P("tIy")))) + P("y0");
			}
		}
		else
		{
			double tshift = 1009;
			double Kyp = Ky;//3390.6+150.0;
			double Kcp = Kc;//1940.6+571.0;
			double dyp = 30.89071766;
			double dcp = 8.924201503;
			double yp = 2274.133385;
			double cp = 1659.959446;
			double p0 = P("p0");
			double tIc = P("tIc");
			double tIy = P("tIy");
			double ryp = dyp*Kyp/yp/(Kyp-yp);
			double tIyp = (dyp*Kyp*tshift+yp*Kyp*Math.log((Kyp-yp)/yp)-yp*yp*Math.log((Kyp-yp)/yp))/dyp/Kyp;
			// find splicing params for consumption function
			double rcp = dcp*Kcp/cp/(Kcp-cp);
			double tIcp = (dcp*Kcp*tshift+cp*Kcp*Math.log((Kcp-cp)/cp)-cp*cp*Math.log((Kcp-cp)/cp))/dcp/Kcp;
			// initialize the input functions c and y
			
			System.out.println("ry: " + ryp + "  tIy: " + tIyp);
			System.out.println("rc: " + rcp + "  tIc: " + tIcp);
			
			//historical fitting parameters
			//double Kch = 1940.6;//1609.8;//1925.1;
			//double rch = 0.018671;//0.023148;//0.019399;;
			//double tIch = 10.8366+985;//-12.7237+985;//-17.5155+980;
			//double c0h = 571;//530;//328.0;

			//double Kyh = 3390.6;//9162.8;//5050.5;//3547.6;//3447.6;
			//double ryh = 0.038934;//0.03637;//0.036148;//0.03719;//0.037202;
			//double tIyh = 10.7177+985;//5.8748+985;//12.4573+985;//11.6055+985;//7.3933+985;
			//double y0h = 150;//85;//80.0;
			
			//XYSeries cseries = new XYSeries("C");
			//XYSeries yseries = new XYSeries("Y");
			//XYSeries pseries = new XYSeries("P");
			
			//XYSeries CPseries = new XYSeries("CP");
			//XYSeries Yseries = new XYSeries("YT");
			for (int t=0; t<tspan; t++)//remember origin of time is 600, so 410 corresponds to end of data, i.e. 2009
			{
				if (t<tshift)
				{
					c[t] = P("Kc") / (1 + Math.exp(-P("rc")*(t-P("tIc")))) + P("c0");
					//p[t] = P("Kp") / (1 + Math.exp(-P("rp")*(t-P("tIp")))) + P("p0");
					y[t] = P("Ky") / (1 + Math.exp(-P("ry")*(t-P("tIy")))) + P("y0");
					//c[t] = P("Kc") / (1 + Math.exp(-P("rc")*(t-P("tIc")))) + c0;
					//p[t] = Kp / (1 + Math.exp(-rp*(t+600-tIp))) + p0;
					//y[t] = P("Ky") / (1 + Math.exp(-P("ry")*(t-P("tIy")))) + y0;
				}
				else
				{
					//c[t] = 8.922810900*(t-410) + 1668.852953; //linear
					//y[t] = 30.89071766*(t-410) + 2304.868870; //linear
					c[t] = Kcp / (1 + Math.exp(-rcp*(t-tIcp)));
					////p[t] = Kp / (1 + Math.exp(-rp*(t+1+tshift-tIp)))+p0;//not spliced
					y[t] = Kyp / (1 + Math.exp(-ryp*(t-tIyp)));
				}
				//p[t] = P("Kp") / (1 + Math.exp(-P("rp")*(t-P("tIp"))));
				p[t] = P("Kp") / (1 + Math.exp(-P("rp")*(t-P("tIp")))) + P("p0");
				//System.out.println(t+" "+c[t]+" "+y[t]);
				
				//cseries.add(t,c[t]/Kc);
				//yseries.add(t,y[t]/Ky);
				//pseries.add(t,p[t]/Kp);
				
				//CPseries.add(t,c[t]*p[t]/T);
				//Yseries.add(t,y[t]);
					
			}
		}
	}
	
	private void runSimulation()
	{
		for (int t=0; t<tspan-1; t++)
		{
			double Rplus = c[t+1]*p[t+1]/y[t+1]-AF[t]-AP[t];
			double Rminus = -Rplus;
			
			F[t+1] = F[t] - P("alpha")*Rplus*Theta(Rplus) + P("beta")*BF[t] - P("zeta")*P("s")*(p[t+1]-p[t]);
			P[t+1] = P[t] - (1-P("alpha"))*Rplus*Theta(Rplus) + P("delta")*BP[t] - (1-P("zeta"))*P("s")*(p[t+1]-p[t]);
			AF[t+1] = AF[t] + P("alpha")*Rplus*Theta(Rplus) - P("gamma")*Rminus*Theta(Rminus) - P("epsilon")*AF[t];
			AP[t+1] = AP[t] + (1-P("alpha"))*Rplus*Theta(Rplus) - (1-P("gamma"))*Rminus*Theta(Rminus) - P("epsilon")*AP[t];
			BF[t+1] = BF[t] + P("gamma")*Rminus*Theta(Rminus) - P("beta")*BF[t] + P("epsilon")*AF[t];
			BP[t+1] = BP[t] + (1-P("gamma"))*Rminus*Theta(Rminus) - P("delta")*BP[t] + P("epsilon")*AP[t];
			U[t+1] = U[t] + P("s")*(p[t+1]-p[t]);
		}
	}
	
	private void runPhaseAnalysis() throws IOException
	{	
		double Kymin = 2400;
		double Kymax = 10000;
		double Kcmin = 1700;
		double Kcmax = 6000;
		
		double ry,tIy,rc,tIc;
		
		double dyp = 30.89071766;
		double dcp = 8.924201503;
		double yp = 2274.133385;
		double cp = 1659.959446;
		double p0 = P("p0");
		//double tIc = P("tIc");
		//double tIy = P("tIy");
		double tIp = P("tIp");
		double rp = P("rp");
		double Kp = P("Kp");
		double s = P("s");
		double alpha = P("alpha");
		double beta = P("beta");
		double gamma = P("gamma");
		double delta = P("delta");
		double epsilon = P("epsilon");
		double zeta = P("zeta");
		
		double tshift = 1009;
		
		double iter = 0;
		
		PrintWriter dataWriter = createFileWriter("phasediagram.dat");
		
		for (double Ky=Kymin; Ky <= Kymax; Ky=Ky+50)
		{
			// find splicing params for yield function
			ry = dyp*Ky/yp/(Ky-yp);
			tIy = (dyp*Ky*tshift+yp*Ky*Math.log((Ky-yp)/yp)-yp*yp*Math.log((Ky-yp)/yp))/dyp/Ky;
			for (double Kc=Kcmin; Kc<=Kcmax; Kc=Kc+25)
			{
				iter++;//useless.
				System.out.println(Ky + "\t" + Kc);
				// find splicing params for consumption function
				rc = dcp*Kc/cp/(Kc-cp);
				tIc = (dcp*Kc*tshift+cp*Kc*Math.log((Kc-cp)/cp)-cp*cp*Math.log((Kc-cp)/cp))/dcp/Kc;
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
				dataWriter.println(getMetrics(Ky,Kc));
				//if (iter==3070)
				//{
				//	outputResults();
				//}
			}
		}
		dataWriter.close();
	}
	
	private String getMetrics(double Ky, double Kc) throws IOException
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
		
		//dataWriter.println(Ky+"\t"+Kc+"\t"+numFTs+"\t"+numPTs+"\t"+numATs+"\t"+tFT+"\t"+tPT+"\t"+tAT+"\t"+FFT+"\t"+PPT+"\t"+AAT+"\t"+tFequib+"\t"+tPequib+"\t"+tAequib+"\t"+Fequib+"\t"+Pequib+"\t"+Aequib+"\t"+numfFTs+"\t"+numfPTs+"\t"+numfATs+"\t"+tfFT+"\t"+tfPT+"\t"+tfAT+"\t"+FfFT+"\t"+PfPT+"\t"+AfAT+"\t"+phaseflag);
		//dataWriter.close();
		return Ky+"\t"+Kc+"\t"+numFTs+"\t"+numPTs+"\t"+numATs+"\t"+tFT+"\t"+tPT+"\t"+tAT+"\t"+FFT+"\t"+PPT+"\t"+AAT+"\t"+tFequib+"\t"+tPequib+"\t"+tAequib+"\t"+Fequib+"\t"+Pequib+"\t"+Aequib+"\t"+numfFTs+"\t"+numfPTs+"\t"+numfATs+"\t"+tfFT+"\t"+tfPT+"\t"+tfAT+"\t"+FfFT+"\t"+PfPT+"\t"+AfAT+"\t"+phaseflag;
	}
	
	private void runSensitivityAnalysis(double maxDiff, int numSteps) throws IOException
	{
		PrintWriter dataWriter = createFileWriter("sensitivity_"+maxDiff+"_"+numSteps+".dat");
		dataWriter.println("Param"+"\t"+"% Change in param"+"\t"+"Fdiff"+"\t"+"Pdiff"+"\t"+"Adiff"+"\t"+"Udiff");
		//Run the baseline
		runSimulation();
		double FBase = F[F.length-1]/T;
		double PBase = P[P.length-1]/T;
		double ABase = (AF[AF.length-1]+AP[AP.length-1])/T;
		double UBase = U[U.length-1]/T;
			
		String key;
		double value;
		Enumeration keys = params.keys();
		double stepSize = maxDiff/numSteps/100.0;
		while(keys.hasMoreElements()) {
			key = (String)keys.nextElement();
			value = (double)params.get(key);
			for (int i=-numSteps; i<=numSteps; i++)
			{
				double newValue = value*(1.0+i*stepSize);
				params.put(key,newValue);
				initializeInputFns(false,0.0,0.0);
				runSimulation();
				double FVar = F[F.length-1]/T;
				double PVar = P[P.length-1]/T;
				double AVar = (AF[AF.length-1]+AP[AP.length-1])/T;
				double UVar = U[U.length-1]/T;
						
				/*double Fdiff = Math.abs(100.0*(FVar-FBase)/FBase);
				double Pdiff = Math.abs(100.0*(PVar-PBase)/PBase);
				double Adiff = Math.abs(100.0*(AVar-ABase)/ABase);
				double Udiff = Math.abs(100.0*(UVar-UBase)/UBase);*/
				
				double Fdiff = FVar-FBase;
				double Pdiff = PVar-PBase;
				double Adiff = AVar-ABase;
				double Udiff = UVar-UBase;
				
				double diff = Math.sqrt((Fdiff*Fdiff+Pdiff*Pdiff+Adiff*Adiff+Udiff*Udiff)/4.0);
				
				dataWriter.println(key+"\t"+(100.0*i*stepSize)+"\t"+Fdiff+"\t"+Pdiff+"\t"+Adiff+"\t"+Udiff);
				
				System.out.println(key+" "+value+" "+newValue+" "+(i*stepSize)+" "+Adiff);
			}
			params.put(key,value);
		}
		dataWriter.close();
	}
	
	private void outputResults(String filename) throws IOException
	{
		PrintWriter dataWriter = createFileWriter(filename);
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
	
	private PrintWriter createFileWriter(String filename) throws IOException
	{
		PrintWriter dataWriter;
		//try
		//{
			// initialize the file writer and dataset
			String pwd = System.getProperty("user.dir");
			dataWriter = new PrintWriter(
					new BufferedWriter(
							new FileWriter(pwd+"/data/"+filename)));
		//}
		//catch (IOException e1)
		//{
		//	System.out.println("Error creating "+filename);
		//	dataWriter = new PrintWriter();
		//	
		//}
		return dataWriter;
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
	
	private double P(String key)
	{
		return (double)params.get(key);
	}
}
