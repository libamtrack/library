package amtrackservice.server;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.GregorianCalendar;
import java.util.HashMap;

import amtrackservice.client.CalculationService;
import amtrackservice.client.Logger;

import com.google.gwt.user.server.rpc.RemoteServiceServlet;



/**
 * The server side implementation of the RPC service.
 */
@SuppressWarnings("serial")
public class CalculationServiceImpl extends RemoteServiceServlet implements
		CalculationService {
	
	private GregorianCalendar cal;

	public long startCalculation(HashMap<String, String> input,
			String wrapperPath, String functionName)
			throws IllegalArgumentException {
		
		System.out.println("Server side: start calculation using " + wrapperPath + " and " + functionName);
		
		String db_ip="localhost";
		String db_name=getServletContext().getInitParameter("dbname");
		String db_user=getServletContext().getInitParameter("dbuser");
		String db_pass=getServletContext().getInitParameter("dbpass");
		
		CalculationDB db = new CalculationDB(db_ip, db_name, db_user, db_pass);
		cal = new GregorianCalendar();

		long time = cal.getTimeInMillis();
		File file = new File(wrapperPath);
		String resultPath = file.getParentFile().getAbsolutePath() + time;

		StringBuffer fileText = new StringBuffer("");
		for (String s : input.keySet()) {
			fileText.append(s + ":" + input.get(s) + "\n");
		}
		try {
			CalculationFile.writeData(new String(fileText), resultPath);
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
		System.out.println("Server side: inserted into DB " + time + " and " + resultPath);
		return time;
	}

	@Override
	public HashMap<String, String> requestResult(long id)
			throws IllegalArgumentException {

		System.out.println("Server side: getting results for id = " + id );
		
		String db_ip="localhost";
		String db_name=getServletContext().getInitParameter("dbname");
		String db_user=getServletContext().getInitParameter("dbuser");
		String db_pass=getServletContext().getInitParameter("dbpass");
		
		HashMap<String, String> result = new HashMap<String, String>();
		CalculationDB db = new CalculationDB(db_ip, db_name, db_user, db_pass);
		String calcPath = db.getCalculation(id).getFilePath();
		try {
			result.putAll(CalculationFile.readData(calcPath));
		} catch (FileNotFoundException e) {
			Logger.error("Results file " + calcPath + " not ready yet");
		} catch (IOException e) {
			e.printStackTrace();
		}

		return result;

	}

	@Override
	public HashMap<Long, String> getCalculations()
			throws IllegalArgumentException {
		
		String db_ip="localhost";
		String db_name=getServletContext().getInitParameter("dbname");
		String db_user=getServletContext().getInitParameter("dbuser");
		String db_pass=getServletContext().getInitParameter("dbpass");
		
		CalculationDB db = new CalculationDB(db_ip, db_name, db_user, db_pass);
		CalculationData[] data = db.getAllCalculations();
		HashMap<Long, String> calculations = new HashMap<Long, String>();
		for (CalculationData d : data) {
			calculations.put(d.getId(), d.getFunctionName());
		}
		
		System.out.println("Server side: get all calculation (number) = " + calculations.size() );
		
		return calculations;
	}
}
