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

	public CalculationDB(String dbIP, String dbName, String username,
			String password) {
		this.url = "jdbc:mysql://" + dbIP + "/";
		this.db = dbName;
		this.driver = "com.mysql.jdbc.Driver";
		this.user = username;
		this.pass = password;
		this.url = this.url + this.db;
		
	}
	public CalculationData getCalculation(long id) {

		String resultPath ="";
		String functionName = "";
		boolean calculated = false;

		String query = "SELECT Path,running,function FROM `"+db+"`.`calculations` where `ID` = ?";
		try {
			Connection connection = getConn(url,user,pass,driver);
			PreparedStatement prep = connection.prepareStatement(query);
			prep.setLong(1, id);
			ResultSet result = prep.executeQuery();

			result.next();
			resultPath = result.getString(1);
			calculated = result.getBoolean(2);
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

		String query = "SELECT * FROM  `amtrack`.`calculations`";

		// prepare for rpc transport
		CalculationData[] calcs = null;

		try {
			Connection connection = getConn(url,user,pass,driver);
			ResultSet result = connection.prepareStatement(query).executeQuery();

			int rsSize = getResultSetSize(result); //size the array
			calcs = new CalculationData[rsSize];

			int i = 0;
			while (result.next()) {
				calcs[i] = new CalculationData(result.getLong(1),result.getString(2),result.getBoolean(3),result.getString(4));
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

		String query = "INSERT INTO `amtrack`.`calculations` ( `ID` , `Path`, `function` ) VALUES ( ?, ?, ? );";
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

		String query = "UPDATE  `amtrack`.`calculations` SET  `running` =  0 WHERE  `calculations`.`ID` = ?;";

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
