const world = new ji.World(0, 9.8 * 100);
const car = new ji.RigidBody({
    pos: ji.V.V(0, 0),
    mass: 100,
    inertia: 100000,
    friction: 0.8,
    shape: new ji.shape.Box(0, 0, 0, 150, 50)
});
const carTop = new ji.RigidBody({
    pos: ji.V.V(0, -50 / 2 - 25),
    mass: 100,
    inertia: 100000,
    friction: 0.8,
    shape: new ji.shape.Convex([
        ji.V.V(-30, -25), ji.V.V(30, -25),
        ji.V.V(50, 25), ji.V.V(-50, 25)
    ])
});
const leftWheel = new ji.RigidBody({
    pos: ji.V.V(-45, 50 / 2),
    mass: 100,
    inertia: 50000,
    friction: 1,
    shape: new ji.shape.Circle(0, 0, 20)
});
const rightWheel = new ji.RigidBody({
    pos: ji.V.V(45, 50 / 2),
    mass: 100,
    inertia: 50000,
    friction: 1,
    shape: new ji.shape.Circle(0, 0, 20)
});
world.add(car, "car");
world.add(carTop, "car");
world.add(leftWheel, "wheels");
world.add(rightWheel, "wheels");
const terrain = {};
const cam = { x: 0, y: 0 };

let keys = {};
function keyPressed() {
    keys[keyCode] = true;
    keys[key.toString().toLowerCase()] = true;
}
function keyReleased() {
    keys[keyCode] = false;
    keys[key.toString().toLowerCase()] = false;
}

let anchor = false, anchorBody;

function generateTerrain(i) {
    return noise(i * 0.01 + 100) * 50 + noise(i * 0.001 + 100) * 500 + noise(i * 0.0001 + 100) * 1000 + height;
}

let pt = Date.now();
function draw() {
    const t = Date.now();
    const dt = constrain(t - pt, 0, 20);
    const cmouseX = mouseX + cam.x - width / 2, cmouseY = mouseY + cam.y - height / 2;

    const terrainInc = 50;
    const start = Math.floor((cam.x - width / 2) / terrainInc) * terrainInc, end = Math.ceil((cam.x + width / 2) / terrainInc) * terrainInc;
    const bottomTerrain = 200 + generateTerrain(0);
    for(let i = start; i < end; i += terrainInc) {
        if(!terrain[i]) {
            const block = new ji.RigidBody({
                pos: ji.V.V(i + terrainInc / 2, 0),
                mass: Infinity, inertia: Infinity,
                friction: 1,
                shape: new ji.shape.Convex([
                    ji.V.V(-terrainInc / 2, bottomTerrain), ji.V.V(terrainInc / 2, bottomTerrain),
                    ji.V.V(terrainInc / 2, bottomTerrain - generateTerrain(i + terrainInc)), ji.V.V(-terrainInc / 2, bottomTerrain - generateTerrain(i))
                ]),
                static: true
            });
            world.add(block);
            terrain[i] = block;
        }
    }
    for(const block in terrain) {
        if(terrain[block].pos.x < start || terrain[block].pos.x > end + terrainInc) {
            world.remove(terrain[block]);
            delete terrain[block];
        }
    }

    background(0, 200, 255);
    const iterations = 6;
    for(let i = 0; i < iterations; i += 1) {
        if(keys[LEFT_ARROW] || keys.a) {
            leftWheel.rvel -= dt * 0.02;
            rightWheel.rvel -= dt * 0.02;
        }
        if(keys[RIGHT_ARROW] || keys.d) {
            leftWheel.rvel += dt * 0.02;
            rightWheel.rvel += dt * 0.02;
        }
        if(keys[UP_ARROW] || keys.w) {
            car.rvel += dt * 0.05;
        }
        if(keys[DOWN_ARROW] || keys.s) {
            car.rvel -= dt * 0.05;
        }
        leftWheel.rvel -= leftWheel.rvel * 0.0005 * dt;
        rightWheel.rvel -= rightWheel.rvel * 0.0005 * dt;
        world.simulate(dt * 0.001 / iterations);
    }
    world.clearConstraints();
    world.addConstraint("union", {
        body1: car,
        body2: carTop,
        pos: ji.V.V(0, -50 / 2 - 25),
        rot: 0
    });
    world.addConstraint("joint", {
        body1: car,
        anchor1: ji.V.V(-45, 50 / 2 + 20),
        body2: leftWheel,
        anchor2: ji.V.V(0, 0)
    // }, "hard", 0.3, 0.05);
    // }, "oscillation", Math.sqrt(0.3 / (4 * Math.PI * Math.PI * 0.05 * 0.002777777778)), 2 * Math.log(2) * 0.05 / (1 - 0.3));
    }, "oscillation", 7.4, 0.1);
    world.addConstraint("joint", {
        body1: car,
        anchor1: ji.V.V(45, 50 / 2 + 20),
        body2: rightWheel,
        anchor2: ji.V.V(0, 0)
    }, "oscillation", 7.4, 0.1);
    noFill(255);
    stroke(0);
    strokeWeight(1);
    if(mouseIsPressed && anchor) {
        // world.addConstraint("anchor", {
        //     body: anchorBody,
        //     anchor,
        //     point: ji.V.V(cmouseX, cmouseY),
        // });
    }
    push();
    translate(width / 2 - cam.x, height / 2 - cam.y);
    for(const body of world.bodies) {
        const shape = body.shape;
        noFill();
        stroke(0);
        if(body === car) fill(255, 0, 0);
        else if(body === carTop) fill(0, 0, 255);
        else if(body !== leftWheel && body !== rightWheel) {
            stroke(0, 200, 50);
            fill(0, 200, 50);
        }
        if(!anchor && body.pointInside(ji.V.V(cmouseX, cmouseY))) {
            if(mouseIsPressed) {
                anchor = ji.V.mdot(body.itransform(), ji.V.V(cmouseX, cmouseY));
                anchorBody = body;
            }
        }
        push();
        translate(body.pos.x, body.pos.y);
        rotate(body.rot);
        if(shape instanceof ji.shape.Box) {
            rect(-shape.size.x / 2, -shape.size.y / 2, shape.size.x, shape.size.y);
        }
        else if(shape instanceof ji.shape.Convex) {
            beginShape();
            for(const vert of shape.verts) {
                vertex(vert.x, vert.y);
            }
            endShape(CLOSE);
        }
        else if(shape instanceof ji.shape.Circle) {
            ellipse(shape.pos.x, shape.pos.y, shape.radius * 2, shape.radius * 2);
            for(let i = 0; i < Math.PI * 2; i += Math.PI * 2 / 2) {
                line(0, 0, shape.pos.x + Math.cos(i) * shape.radius, shape.pos.y + Math.sin(i) * shape.radius);
            }
        }
        pop();
    }
    pop();
    cam.x += (car.pos.x - cam.x + width / 3) * 0.2;
    cam.y += (car.pos.y - cam.y) * 0.2;
    if(!mouseIsPressed) anchor = false;
    pt = t;
}