int numBodies;
float[][][] positions; // positions[frame][body][x,y,z]
int totalFrames;
int currentFrame = 0;
float dt = 0.016;
float camDistance = 60;  // distance caméra
float offsetX = 0;
float offsetY = 0;
float lastMouseX, lastMouseY;
boolean dragging = false;


// centrage + échelle
float cx = 0, cy = 0, cz = 0;
float scaleFactor;

float zoom = 1.0; // <1 pour zoom arrière

void mouseWheel(processing.event.MouseEvent event) {
  if(zoom+event.getCount()*0.05 > 1) return;
  zoom += event.getCount()*0.05;
}
void mousePressed() {
  if (mouseButton == LEFT) {
    dragging = true;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
  }
}

void mouseReleased() {
  if (mouseButton == LEFT) {
    dragging = false;
  }
}

void mouseDragged() {
  if (dragging) {
    float dx = mouseX - lastMouseX;
    float dy = mouseY - lastMouseY;
    offsetX += dx;
    offsetY += dy;
    lastMouseX = mouseX;
    lastMouseY = mouseY;
  }
}



void setup() {
  pixelDensity(2);
  size(800, 600, P3D);
  ortho();
  String[] lines = loadStrings("../results/Barnes_Hut_parallelstars_data");
  numBodies = int(lines[0]);
  totalFrames = lines.length - 1;
  positions = new float[totalFrames][numBodies][3];

  // --- Lecture du CSV ---
  for (int f = 0; f < totalFrames; f++) {
    String[] values = split(lines[f+1], ',');
    for (int b = 0; b < numBodies; b++) {
      int idx = b * 3;
      positions[f][b][0] = float(values[idx]);
      positions[f][b][1] = float(values[idx + 1]);
      positions[f][b][2] = float(values[idx + 2]);
    }
  }

  // --- Calcul du centre global ---
  for (int f = 0; f < totalFrames; f++) {
    for (int b = 0; b < numBodies; b++) {
      cx += positions[f][b][0];
      cy += positions[f][b][1];
      cz += positions[f][b][2];
    }
  }

  // --- Calcul du centre global (frame 0 uniquement) ---
  cx = 0;
  cy = 0;
  cz = 0;
  for (int b = 0; b < numBodies; b++) {
    cx += positions[0][b][0];
    cy += positions[0][b][1];
    cz += positions[0][b][2];
  }
  cx /= numBodies;
  cy /= numBodies;
  cz /= numBodies;
  
  // --- Calcul de l'échelle automatique (frame 0, rayon cylindrique) ---
  float maxR = 0;
  for (int b = 0; b < numBodies; b++) {
    float dx = positions[0][b][0] - cx;
    float dy = positions[0][b][1] - cy;
    float R = sqrt(dx*dx + dy*dy);  // rayon dans le plan XY
    maxR = max(maxR, R);
  }
  
  // scaleFactor : plus petit pour zoomer, plus grand pour dézoomer
scaleFactor = maxR / (min(width, height) * 0.45);
}

void draw() {
  background(0);

  // Projection orthographique
  ortho();

  // Centrer l'origine
  translate(width/2, height/2);

  stroke(0, 150, 255, 90);
  strokeWeight(2);
  
  // Tous les points dans un seul beginShape
  beginShape(POINTS);
  for (int b = 0; b < numBodies; b++) {
    float x = (positions[currentFrame][b][0] - cx) / scaleFactor;
    float y = (positions[currentFrame][b][1] - cy) / scaleFactor;
    float z = (positions[currentFrame][b][2] - cz) / scaleFactor;
    
    vertex((x-offsetX*0.05)*zoom, (y-offsetY*0.05)*zoom, z*zoom);
  }
  endShape();

  // Texte 2D
  hint(DISABLE_DEPTH_TEST);
  fill(255);
  textAlign(RIGHT, TOP);
  textSize(24);
  text("dt = " + dt, width - 10, 10);
  hint(ENABLE_DEPTH_TEST);

  currentFrame = (currentFrame + 1) % totalFrames;
  dt += 0.01;
}
