/* -------------------
* ScatterPlot.java
* -------------------
* This class handles the creation and rendering of scatter plots given an XYDataset object and
* frame width and frame height.  It is based on LineChartDemo2.java from the JFreeChart
* Developer Guide
* (C) Copyright 2002-2005, by Object Refinery Limited.
*
*/

//package fpamodel;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.geom.Point2D;
import javax.swing.JPanel;
import org.jfree.chart.*;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.*;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;

// Referenced classes of package demo:
//			SampleXYDataset2

public class ScatterPlot extends ApplicationFrame
{
	private XYDataset dataset;

	static class MyChartMouseListener
		implements ChartMouseListener
	{

		ChartPanel panel;

		public void chartMouseClicked(ChartMouseEvent chartmouseevent)
		{
			int i = chartmouseevent.getTrigger().getX();
			int j = chartmouseevent.getTrigger().getY();
			Point2D point2d = panel.translateScreenToJava2D(new Point(i, j));
			XYPlot xyplot = (XYPlot)panel.getChart().getPlot();
			ChartRenderingInfo chartrenderinginfo = panel.getChartRenderingInfo();
			java.awt.geom.Rectangle2D rectangle2d = chartrenderinginfo.getPlotInfo().getDataArea();
			double d = xyplot.getDomainAxis().java2DToValue(point2d.getX(), rectangle2d, xyplot.getDomainAxisEdge());
			double d1 = xyplot.getRangeAxis().java2DToValue(point2d.getY(), rectangle2d, xyplot.getRangeAxisEdge());
			ValueAxis valueaxis = xyplot.getDomainAxis();
			ValueAxis valueaxis1 = xyplot.getRangeAxis();
			double d2 = valueaxis.valueToJava2D(d, rectangle2d, xyplot.getDomainAxisEdge());
			double d3 = valueaxis1.valueToJava2D(d1, rectangle2d, xyplot.getRangeAxisEdge());
			Point point = panel.translateJava2DToScreen(new java.awt.geom.Point2D.Double(d2, d3));
			System.out.println("Mouse coordinates are (" + i + ", " + j + "), in data space = (" + d + ", " + d1 + ").");
			System.out.println("--> (" + point.getX() + ", " + point.getY() + ")");
		}

		public void chartMouseMoved(ChartMouseEvent chartmouseevent)
		{
		}

		public MyChartMouseListener(ChartPanel chartpanel)
		{
			panel = chartpanel;
		}
	}


	/*public ScatterPlot(String s)
	{
		super(s);
		JPanel jpanel = createDemoPanel();
		jpanel.setPreferredSize(new Dimension(500, 270));
		setContentPane(jpanel);
	}*/
	public ScatterPlot(String title, XYDataset data, int chartWidth, int chartHeight) 
	{
		super(title);
		dataset = data;//createDataset();
		JPanel jpanel = createPanel();
		//JFreeChart chart = createChart(dataset);
		//ChartPanel chartPanel = new ChartPanel(chart);
		jpanel.setPreferredSize(new Dimension(chartWidth, chartHeight));
		setContentPane(jpanel);
	}

	private JFreeChart createChart(XYDataset xydataset)
	{
		JFreeChart jfreechart = ChartFactory.createScatterPlot(null, "X", "Y", xydataset, PlotOrientation.VERTICAL, true, true, false);
		XYPlot xyplot = (XYPlot)jfreechart.getPlot();
		xyplot.setDomainCrosshairVisible(true);
		xyplot.setDomainCrosshairLockedOnData(true);
		xyplot.setRangeCrosshairVisible(true);
		xyplot.setRangeCrosshairLockedOnData(true);
		xyplot.setDomainZeroBaselineVisible(true);
		xyplot.setRangeZeroBaselineVisible(true);
		xyplot.setDomainPannable(true);
		xyplot.setRangePannable(true);
		xyplot.setBackgroundPaint(Color.white);
		//xyplot.setAxisOffset(new RectangleInsets(5.0, 5.0, 5.0, 5.0));
		xyplot.setDomainGridlinePaint(Color.lightGray);
		xyplot.setRangeGridlinePaint(Color.lightGray);
		NumberAxis numberaxis = (NumberAxis)xyplot.getDomainAxis();
		numberaxis.setAutoRangeIncludesZero(false);
		return jfreechart;
	}

	public JPanel createPanel()
	{
		JFreeChart jfreechart = createChart(dataset);
		ChartPanel chartpanel = new ChartPanel(jfreechart);
		chartpanel.setMouseWheelEnabled(true);
		chartpanel.addChartMouseListener(new MyChartMouseListener(chartpanel));
		return chartpanel;
	}
}
