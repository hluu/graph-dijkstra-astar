package roadgraph;

import geography.GeographicPoint;

public class StreetSegment {
	private GeographicPoint point1;
	private GeographicPoint point2;
	
	private String roadName;
	private String roadType;
	
	// Length in km
	private double length;
	
	public StreetSegment(GeographicPoint pt1, GeographicPoint pt2, String roadName,
			String roadType, double length) {
		this.point1 = pt1;
		this.point2 = pt2;
		this.roadName = roadName;
		this.roadType = roadType;
		this.length = length;
	}

	public GeographicPoint getPoint1() {
		return point1;
	}

	public GeographicPoint getPoint2() {
		return point2;
	}

	public String getRoadName() {
		return roadName;
	}

	public String getRoadType() {
		return roadType;
	}

	public double getLength() {
		return length;
	}
	
	
}
