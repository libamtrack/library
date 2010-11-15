package amtrackservice.server;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.PreparedStatement;

public class CalculationDB extends AbstractDBConnection {
	
	private String url = null;
	private String db = null;
	private String driver = null;
	private String user = null;
	private String pass = null;

	public CalculationDB(String dbIP, String dbName, String username, String password) {
		this.url = "jdbc:sqlite:";
		this.db = dbName;
		this.driver = "org.sqlite.JDBC";
		this.user = username;
		this.pass = password;
		this.url = this.url + this.db;
		
	}
	public CalculationData getCalculation(long id) {

		String resultPath ="";
		String functionName = "";
		boolean calculated = false;

		String query = "SELECT Path,running,function FROM calculations where ID = ?";
		try {
			Connection connection = getConn(url,user,pass,driver);
			PreparedStatement prep = connection.prepareStatement(query);
			prep.setLong(1, id);
			ResultSet result = prep.executeQuery();
			result.next();
			resultPath = result.getString(1);
			if( result.getLong(2) == 0 )
				calculated = false;
			else
				calculated = true;
			functionName = result.getString(3);

			result.close();
			connection.close();

		} catch (Exception e) {
			System.err.println("Statement Error: " + query);
			e.printStackTrace();
		}

		return new CalculationData(id, resultPath, calculated,functionName);
	}
	
	public CalculationData[] getAllCalculations() {

		String query = "SELECT * FROM calculations";

		// prepare for rpc transport
		CalculationData[] calcs = null;

		try {
			Connection connection = getConn(url,user,pass,driver);
			ResultSet result = connection.prepareStatement(query).executeQuery();

			int rsSize = getResultSetSize(result); //size the array
			calcs = new CalculationData[rsSize];

			int i = 0;
			while (result.next()) {
				boolean calculated;
				if( result.getLong(3) == 0 ){
					calculated = false;
				} else {
					calculated = true;
				}
				calcs[i] = new CalculationData(result.getLong(1),result.getString(2),calculated,result.getString(4));
				i++;
			}

			// clean up
			result.close();
			connection.close();

		} catch(Exception e) {
			System.err.println("Statement Error: " + query);
			e.printStackTrace();
		}

		// return the array
		return calcs;
	}
	
	public void insertCalculation(long id, String path, String name) {

		String query = "INSERT INTO calculations ( ID , Path, function ) VALUES ( ?, ?, ? );";
		try {			
			Connection connection = getConn(url,user,pass,driver);
			PreparedStatement prep = connection.prepareStatement(query);
			prep.setLong(1,id);
			prep.setString(2, path);
			prep.setString(3,name);
			prep.executeUpdate();
			connection.close();

		} catch(Exception e) {
			System.err.println("Statement Error: " + query);
			e.printStackTrace();
		}
	}
	
	public void finishCalculation(long id) {

		String query = "UPDATE calculations SET running =  0 WHERE  ID = ?;";

		try {
			Connection connection = getConn(url,user,pass,driver);
			PreparedStatement prep = connection.prepareStatement(query);
			prep.setLong(1, id);
			prep.executeUpdate();
			connection.close();

		} catch(Exception e) {
			System.err.println("Statement Error: " + query);
			e.printStackTrace();
		}
	}
	
}
