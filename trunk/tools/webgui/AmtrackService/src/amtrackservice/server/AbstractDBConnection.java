package amtrackservice.server;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.ResultSet;
import java.sql.SQLException;

public abstract class AbstractDBConnection {

	protected Connection getConn(String url, String user, String pass,
			String driver) {

		Connection conn = null;

		try {

			Class.forName(driver).newInstance();
			conn = DriverManager.getConnection(url, user, pass);

		} catch (Exception e) {
			System.err.println("Connection Error: ");
			e.printStackTrace();
		}

		if (conn == null) {
			System.out.println("Can't connect to MyDatabase");
		}

		return conn;
	}

	protected static int getResultSetSize(ResultSet resultSet) {
		int size = -1;

		try {
			resultSet.last();
			size = resultSet.getRow();
			resultSet.beforeFirst();
		} catch (SQLException e) {
			return size;
		}

		return size;
	}

}