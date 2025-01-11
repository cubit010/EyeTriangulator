#include <iostream>
#include <cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <regex>
#include <vector>
#include <algorithm>
#include <windows.h>
#include <limits>  // For std::numeric_limits

// Constants
const double PI = 3.14159265358979323846;

struct Point {
    double x;
    double z;
};

// Convert degrees to radians
double toRadians(double degrees) {
    return degrees * PI / 180.0;
}

// Normalize angle to 0-360 degrees
double normalizeAngle(double angle) {
    angle = 270.00 - angle;
    return angle;
}

// Normalize angle tester
double normalizeAngletest(double angle) {
    double angleTest = 270.00 - angle;
    while (angleTest < 0 || angleTest > 360) {
        if (angleTest < 0) {
            angleTest += 360;
        }
        else if (angleTest > 360) {
            angleTest -= 360;
        }
    }
    return angleTest;
}

// Calculate the intersection point of two lines defined by points and angles
bool triangulate(const Point& p1, double angle1, const Point& p2, double angle2, Point& intersection) {
    // Normalize angles
    angle1 = normalizeAngle(angle1);
    angle2 = normalizeAngle(angle2);

    // Convert angles to radians
    double theta1 = toRadians(angle1);
    double theta2 = toRadians(angle2);

    // Line equations: y = mx + b
    // For line 1
    double m1 = std::tan(theta1);
    double b1 = p1.z - m1 * p1.x;

    // For line 2
    double m2 = std::tan(theta2);
    double b2 = p2.z - m2 * p2.x;

    // Check for parallel lines
    if (std::abs(m1 - m2) < 1e-9) {
        return false; // No intersection (lines are parallel)
    }

    // Calculate intersection
    intersection.x = (b2 - b1) / (m1 - m2);
    intersection.z = m1 * intersection.x + b1;
    return true;
}

// Calculate the distance between two points
double calculateDistance(const Point& p1, const Point& p2) {
    return std::sqrt(std::pow(p2.x - p1.x, 2) + std::pow(p2.z - p1.z, 2));
}

// Find the furthest vertex from the centroid
void findFurthestVertex(const std::vector<Point>& points, const Point& centroid, Point& furthestPoint, double& maxDistance) {
    maxDistance = 0;
    for (const auto& point : points) {
        double distance = calculateDistance(point, centroid);
        if (distance > maxDistance) {
            maxDistance = distance;
            furthestPoint = point;
        }
    }
}

// Calculate the centroid (center) of a polygon
Point calculatePolygonCentroid(const std::vector<Point>& points) {
    double sumX = 0, sumZ = 0;

    // Sum all the x and z coordinates
    for (const auto& p : points) {
        sumX += p.x;
        sumZ += p.z;
    }

    // Average the coordinates to find the center
    Point center = { sumX / points.size(), sumZ / points.size() };
    return center;
}

// Sort points around a centroid for consistent ordering
bool compareAngle(const Point& a, const Point& b, const Point& centroid) {
    double angleA = std::atan2(a.z - centroid.z, a.x - centroid.x);
    double angleB = std::atan2(b.z - centroid.z, b.x - centroid.x);
    return angleA < angleB;
}

// Sort vertices based on their angle relative to the centroid
std::vector<Point> sortVertices(const std::vector<Point>& points) {
    // Find the centroid
    Point centroid = calculatePolygonCentroid(points);

    // Copy points and sort by angle relative to the centroid
    std::vector<Point> sortedPoints = points;
    std::sort(sortedPoints.begin(), sortedPoints.end(),
        [&centroid](const Point& a, const Point& b) { return compareAngle(a, b, centroid); });

    return sortedPoints;
}

// Calculate the area of a polygon given its vertices using the shoelace formula
double calculatePolygonArea(const std::vector<Point>& points) {
    if (points.size() < 3) return 0.0; // Not a polygon

    // Ensure vertices are ordered
    std::vector<Point> sortedPoints = sortVertices(points);
    std::cout << "Verticies of area of possibilities:\n";
    std::cout << sortedPoints[0].x << ", " << -sortedPoints[0].z << "\n";
	std::cout << sortedPoints[1].x << ", " << -sortedPoints[1].z << "\n";
	std::cout << sortedPoints[2].x << ", " << -sortedPoints[2].z << "\n";
	std::cout << sortedPoints[3].x << ", " << -sortedPoints[3].z << "\n";
    // Compute area using the Shoelace formula
    double area = 0;
    int n = sortedPoints.size();
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n; // Next vertex (looping back to 0 when i is the last vertex)
        area += sortedPoints[i].x * sortedPoints[j].z - sortedPoints[j].x * sortedPoints[i].z;
    }
    return std::abs(area) / 2.0;
}

// Handle input parsing
void parseInput(const std::string& input, Point& point, double& angle) {
    std::regex pattern(R"([-\d\.]+)");
    std::sregex_iterator it(input.begin(), input.end(), pattern), end;

    if (it != end) point.x = std::stod((it++)->str());
    if (it != end) it++; // Skip y-coordinate
    if (it != end) point.z = -std::stod((it++)->str());
    if (it != end) angle = std::stod((it++)->str());
}

// Get clipboard text
std::string getClipboardText() {
    if (!OpenClipboard(nullptr)) return "";
    HANDLE hData = GetClipboardData(CF_TEXT);
    if (hData == nullptr) {
        CloseClipboard();
        return "";
    }
    char* pszText = static_cast<char*>(GlobalLock(hData));
    if (pszText == nullptr) {
        CloseClipboard();
        return "";
    }
    std::string text(pszText);
    GlobalUnlock(hData);
    CloseClipboard();
    return text;
}

// Key listener function
void keyListener(Point& p1, double& angle1, Point& p2, double& angle2) {
    bool isFirstInput = true;

    while (true) {
        if ((GetAsyncKeyState(VK_F3) & 0x8000) && (GetAsyncKeyState(0x43) & 0x8000)) { // Check if F3 and C are pressed simultaneously
            std::string clipboardText = getClipboardText();
            if (clipboardText.empty()) {
                std::cout << "Clipboard is empty or inaccessible.\n";
                continue;
            }

            if (isFirstInput) {
                parseInput(clipboardText, p1, angle1);
                std::cout << "First coordinate: (" << p1.x << ", " << -p1.z << ")\n";
                std::cout << "First angle: " << std::fixed << std::setprecision(2) << normalizeAngletest(angle1) << " degrees\n\n";
                isFirstInput = false;
            }
            else {
                parseInput(clipboardText, p2, angle2);
                std::cout << "Second coordinate: (" << p2.x << ", " << -p2.z << ")\n";
                std::cout << "Second angle: " << std::fixed << std::setprecision(2) << normalizeAngletest(angle2) << " degrees\n\n";
                break; // Exit after second input
            }

            Sleep(300); // Prevent multiple detections
        }
        Sleep(50); // Reduce CPU usage
    }
}

// Function to calculate confidence level
int detConfidence(double area, double radius) {
    double confidence1 = ((2*radius * 2*radius) / area) * 100;
    return confidence1;
}

int confidence95(double area) {
    double confidence = 0;
    int i = 0;
    while (confidence < 95) {
        confidence = detConfidence(area, i);
        i++;
    }
    return i;
}



int main() {
    char continueRunning;
    do {
        Point p1, p2, intersection;
        double angle1, angle2;

        std::cout << "Press F3 + C to input coordinates and angles from clipboard.\n";
        keyListener(p1, angle1, p2, angle2);

        std::vector<Point> intersections;
        double deviation = 0.3;
        double angles1[] = { angle1 - deviation, angle1 + deviation };
        double angles2[] = { angle2 - deviation, angle2 + deviation };

        for (double a1 : angles1) {
            for (double a2 : angles2) {
                Point temp;
                if (triangulate(p1, a1, p2, a2, temp)) {
                    intersections.push_back(temp);
                }
            }
        }

        if (triangulate(p1, angle1, p2, angle2, intersection)) {
            std::cout << "Overworld coordinates are: (" << std::setprecision(0) << intersection.x << ", " << std::setprecision(0) << -intersection.z << ")\n";
            std::cout << "Nether coordinates are: (" << std::setprecision(0) << intersection.x / 8 << ", " << std::setprecision(0) << -intersection.z / 8 << ")\n\n";
        }

        if (intersections.size() == 4) {
            double quadrilateralArea = calculatePolygonArea(intersections);
            Point quadrilateralCenter = calculatePolygonCentroid(intersections);

            Point furthestPoint;
            double maxDistance;
            findFurthestVertex(intersections, quadrilateralCenter, furthestPoint, maxDistance);
			double maxX = intersections[0].x;
			double minX = intersections[0].x;
			double maxZ = intersections[0].z;
			double minZ = intersections[0].z;
            std::cout << "Quadrilateral area formed by two lines' deviations: " << std::fixed << std::setprecision(2) << quadrilateralArea << " square blocks\n";
            std::cout << "Center of the quadrilateral: (" << quadrilateralCenter.x << ", " << -quadrilateralCenter.z << ")\n";
            std::cout << "Radius to search with 95% certainty: " << confidence95(quadrilateralArea) << "\n";
            std::cout << "Max radius to search from center: " << std::fixed << std::setprecision(2) << maxDistance << " blocks\n";
            // Bounding box finding
			for (int i = 0; i < 4; i++) {
				if (intersections[i].x > maxX) {
					maxX = intersections[i].x;
				}
				if (intersections[i].x < minX) {
					minX = intersections[i].x;
				}
				if (intersections[i].z > maxZ) {
					maxZ = intersections[i].z;
				}
				if (intersections[i].z < minZ) {
					minZ = intersections[i].z;
				}
                
			}
			std::cout << "Bounding box: between " << minX << " and " << maxX << " on the X axis, between " << -minZ << " and " << -maxZ << "on the Z axis\n\n";
        }
        else {
            std::cout << "Insufficient intersections to form a polygon.\n";
        }

        std::cout << "Would you like to run again? (y/n): ";
        std::cin >> continueRunning;
        std::cout << std::endl;
        std::cin.ignore(); // Clear newline character
    } while (continueRunning == 'y' || continueRunning == 'Y');

    return 0;
}
