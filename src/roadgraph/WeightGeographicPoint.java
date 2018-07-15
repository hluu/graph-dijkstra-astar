package roadgraph;

import geography.GeographicPoint;

/**
 * This class is used to store both WeightGeographicPoint and weight so we can
 * order them in the priority queue based on the weight.
 * 
 * @author hluu
 *
 */
public class WeightGeographicPoint {
	private GeographicPoint geographicPoint;
	private double weight;
	private double predictedDist;
	
	public static WeightGeographicPoint create(GeographicPoint geographicPoint, double weight) {
		return new WeightGeographicPoint(geographicPoint, weight, 0.0);
	}
	
	public static WeightGeographicPoint create(GeographicPoint geographicPoint, double weight,
			double preditedDistance) {
		return new WeightGeographicPoint(geographicPoint, weight, preditedDistance);
	}
	
	public WeightGeographicPoint(GeographicPoint gp, double weight, double predictedDistance) {
		this.geographicPoint = gp;
		this.weight = weight;
		this.predictedDist = predictedDistance;
	}

	public GeographicPoint getGeographicPoint() {
		return geographicPoint;
	}

	public double getWeight() {
		return weight;
	}
	
	public double getPredictedDist() {
		return predictedDist;
	}
	
	public double getCombinedWeight() {
		return getWeight() + getPredictedDist();
	}
	
	public String toString() {
		return geographicPoint.toString() + " weight: " + weight + " predictedDist: " +  predictedDist;
	}
	
}
