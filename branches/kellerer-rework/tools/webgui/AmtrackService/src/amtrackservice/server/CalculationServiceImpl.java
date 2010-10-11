package amtrackservice.server;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.GregorianCalendar;
import java.util.HashMap;

import amtrackservice.client.CalculationService;
import com.google.gwt.user.server.rpc.RemoteServiceServlet;



/**
 * The server side implementation of the RPC service.
 */
@SuppressWarnings("serial")
public class CalculationServiceImpl extends RemoteServiceServlet implements
		CalculationService {

	private static final String DB_IP="localhost";
	private static final String DB_NAME="";
	private static final String DB_USER="";
	private static final String DB_PASS="";
	
	private GregorianCalendar cal;

	public long startCalculation(HashMap<String, String> input,
			String wrapperPath, String functionName)
			throws IllegalArgumentException {
		CalculationDB db = new CalculationDB(DB_IP, DB_NAME, DB_USER, DB_PASS);
		cal = new GregorianCalendar();

		long time = cal.getTimeInMillis();
		File file = new File(wrapperPath);
		String resultPath = file.getParentFile().getAbsolutePath() + time;

		StringBuffer fileText = new StringBuffer("");
		for (String s : input.keySet()) {
			fileText.append(s + ":" + input.get(s) + "\n");
		}
		try {
			CalculationFile.dateiSchreiben(new String(fileText), resultPath);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		ProcessBuilder pb = new ProcessBuilder(wrapperPath, resultPath);
		try {
			pb.start();
			db.insertCalculation(time, resultPath, functionName);
		} catch (IOException e1) {
			e1.printStackTrace();
		}

		return time;
	}

	@Override
	public HashMap<String, String> requestResult(long id)
			throws IllegalArgumentException {
		HashMap<String, String> result = new HashMap<String, String>();
		CalculationDB db = new CalculationDB(DB_IP, DB_NAME, DB_USER, DB_PASS);
		String calcPath = db.getCalculation(id).getFilePath();
		try {
			result.putAll(CalculationFile.dateiAuslesen(calcPath));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;

	}

	@Override
	public HashMap<Long, String> getCalculations()
			throws IllegalArgumentException {
		CalculationDB db = new CalculationDB(DB_IP, DB_NAME, DB_USER, DB_PASS);
		CalculationData[] data = db.getAllCalculations();
		HashMap<Long, String> calculations = new HashMap<Long, String>();
		for (CalculationData d : data) {
			calculations.put(d.getId(), d.getFunctionName());
		}
		return calculations;
	}
}
