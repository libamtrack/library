package amtrackservice.client;

import java.util.HashMap;

import amtrackservice.client.gui.elements.AmWidget;
import amtrackservice.shared.Config;

import com.google.gwt.core.client.GWT;
import com.google.gwt.user.client.rpc.AsyncCallback;

/**
 * implements methods to contact the server.
 * @author Christoph Kolb, 2010 Hochschule Heilbronn
 *
 */
public class CalculationControl {

	/**
	 * Create a remote service proxy to talk to the server-side service.
	 */
	private final CalculationServiceAsync calculationService = GWT.create(CalculationService.class);

	private Config config = null;
	
	private HashMap<Long,String> runningOnServer = new HashMap<Long, String>();

	private static CalculationControl instance = null;

	private CalculationControl() {
		// TODO Auto-generated constructor stub
	}

	public static CalculationControl getInstance() {
		if (instance == null) {
			instance = new CalculationControl();
		}
		return instance;
	}

	public void setConfig(Config config) {
		this.config = config;
	}

	public void getCalculationResult(final Calculation calculation) {
		Logger.info("Requesting result for calculation " + calculation.getName() + "(" + calculation.getTimeID() + ")");
		Logger.init(calculation);
		calculationService.requestResult(calculation.getTimeID(),
				new AsyncCallback<HashMap<String, String>>() {

					@Override
					public void onSuccess(HashMap<String, String> result) {						
						Logger.info("Reqest successful");
						if( result.get("exitcode") != null){
							String exitcode = result.get("exitcode").trim();
							if( exitcode == ""){
								Logger.info("Process not finished yet");
								calculation.setCalculationResult(new HashMap<String, String>());
							}else{
								Logger.info("Standard output : <BR>" + result.get("stdout").replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
								Logger.info("Standard error : <BR>" + result.get("stderr").replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
								if( exitcode.equals("0")){
									Logger.info("Process finished, success" );
									Logger.info("Output file content : <BR>" + result.get("content").replaceAll("(\r\n|\r|\n|\n\r)", "<br>"));
									calculation.setCalculationResult(result);
								} else {
									Logger.info("Process finished, failure, exit code = _" + exitcode + "_");
									calculation.setCalculationResult(new HashMap<String, String>());
									calculation.showStatus();
								}
							}
						} else {
							Logger.info("Process not finished yet");
							calculation.setCalculationResult(new HashMap<String, String>());
						}
					}

					@Override
					public void onFailure(Throwable caught) {
						Logger.log(caught.getMessage());
					}
				});

	}

	public void calculate(final Calculation calculation) {
		Logger.info("Starting calculation for " + calculation.getName());
		
		HashMap<String, String> calcData = new HashMap<String, String>();
		for (AmWidget widget : calculation.getInputWidgets()) {
			calcData.put(widget.getDataLink(), widget.getValue());
		}

		calculationService.startCalculation(calcData, config.getWrapperPath()
				+ calculation.getWrapper(),calculation.getName(), new AsyncCallback<Long>() {

			@Override
			public void onSuccess(Long result) {
				calculation.setTimeID(result);
				Logger.info(calculation.getName() + ": calculation request sent. Id set to: " + result);
			}

			@Override
			public void onFailure(Throwable caught) {
				Logger.log(caught.getMessage());
			}
		});

	}

	public void getCalculation(final long id, final AmtrackService amtrackService) {
		calculationService.getCalculations(new AsyncCallback<HashMap<Long,String>>() {

			@Override
			public void onFailure(Throwable caught) {
				Logger.log(caught.getMessage());
				
			}

			@Override
			public void onSuccess(HashMap<Long, String> result) {
				result.get(id);
			}
		});	
	}
	
	private void setRunning(HashMap<Long, String> running){
		this.runningOnServer = running;
	}
	
	public HashMap<Long, String> getRunning(){
		return runningOnServer;
	}

}
