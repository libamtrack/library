package amtrackservice.client.gui.elements;

import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.JsArray;
import com.google.gwt.dom.client.Element;
import com.google.gwt.visualization.client.AbstractDataTable;
import com.google.gwt.visualization.client.AbstractDrawOptions;
import com.google.gwt.visualization.client.Color;
import com.google.gwt.visualization.client.CommonOptions;
import com.google.gwt.visualization.client.Selectable;
import com.google.gwt.visualization.client.Selection;
import com.google.gwt.visualization.client.events.Handler;
import com.google.gwt.visualization.client.events.OnMouseOutHandler;
import com.google.gwt.visualization.client.events.OnMouseOverHandler;
import com.google.gwt.visualization.client.events.ReadyHandler;
import com.google.gwt.visualization.client.events.SelectHandler;
import com.google.gwt.visualization.client.visualizations.Visualization;

public class ScatterChartLogAxis extends Visualization<ScatterChartLogAxis.Options> implements
Selectable {

	public static class AxisOptions extends AbstractDrawOptions {
		public static AxisOptions create() {
			return JavaScriptObject.createObject().cast();
		}

		protected AxisOptions() {
		}

		public final native void setBaseline(double baseline) /*-{
		    this.baseline = baseline;
		  }-*/;

		public final native void setBaselineColor(String baselineColor) /*-{
		    this.baselineColor = baselineColor;
		  }-*/;

		public final native void setIsLogScale(boolean isLogScale) /*-{
		    this.logScale = isLogScale;
		  }-*/;

		public final native void setMaxValue(double max) /*-{
		    this.maxValue = max;
		  }-*/;

		public final native void setMinValue(double min) /*-{
		    this.minValue = min;
		  }-*/;

		public final native void setTitle(String title) /*-{
		    this.title = title;
		  }-*/;

	}


	/**
	 * Options for drawing the chart.
	 */
	public static class Options extends CommonOptions {
		public static Options create() {
			return JavaScriptObject.createObject().cast();
		}

		protected Options() {
		}

		public final native void setLineSize(int size) /*-{
      this.lineSize = size;
    }-*/;

		public final native void setPointSize(int size) /*-{
      this.pointSize = size;
    }-*/;

		public final native void setAxisBackgroundColor(Color color) /*-{
    this.axisBackgroundColor = color;
  }-*/;

		public final native void setAxisBackgroundColor(String color) /*-{
    this.axisBackgroundColor = color;
  }-*/;

		public final native void setAxisColor(Color color) /*-{
    this.axisColor = color;
  }-*/;

		public final native void setAxisColor(String color) /*-{
    this.axisColor = color;
  }-*/;

		public final native void setAxisFontSize(double size) /*-{
    this.axisFontSize = size;
  }-*/;

		public final native void setMax(double max) /*-{
    this.max = max;
  }-*/;

		public final native void setMin(double min) /*-{
    this.min = min;
  }-*/;

		public final native void setReverseAxis(boolean reverseAxis) /*-{
    this.reverseAxis = reverseAxis;
  }-*/;

		public final native void setShowCategories(boolean showCategories) /*-{
            this.showCategories = showCategories;
        }-*/;

		public final native void setTitleX(String title) /*-{
            this.titleX = title;
        }-*/;

		public final native void setTitleY(String title) /*-{
            this.titleY = title;
        }-*/;

		public final native void setHAxisOptions(AxisOptions options) /*-{
		    this.hAxis = options;
        }-*/;

		public final native void setVAxisOptions(AxisOptions options) /*-{
	        this.vAxis = options;
        }-*/;

	}

	public static final String PACKAGE = "scatterchartlog";

	public ScatterChartLogAxis() {
		super();
	}

	public ScatterChartLogAxis(AbstractDataTable data, Options options) {
		super(data, options);
	}

	public final void addOnMouseOutHandler(OnMouseOutHandler handler) {
		Handler.addHandler(this, "onmouseout", handler);
	}

	public final void addOnMouseOverHandler(OnMouseOverHandler handler) {
		Handler.addHandler(this, "onmouseover", handler);
	}

	public final void addReadyHandler(ReadyHandler handler) {
		Handler.addHandler(this, "ready", handler);
	}

	public final void addSelectHandler(SelectHandler handler) {
		Selection.addSelectHandler(this, handler);
	}

	public final JsArray<Selection> getSelections() {
		return Selection.getSelections(this);
	}

	public final void setSelections(JsArray<Selection> sel) {
		Selection.setSelections(this, sel);
	}

	@Override
	protected native JavaScriptObject createJso(Element parent) /*-{
    return new $wnd.google.visualization.ScatterChart(parent);
  }-*/;  

}
