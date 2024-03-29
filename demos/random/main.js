const world = new ji.World(0, 9.8 * 100);
const ropeLength = 20;

world.add(new ji.RigidBody({
    pos: ji.V.V(width / 2, 0),
    shape: new ji.shape.Box(0, 0, 0, width, 200),
    mass: Infinity,
    inertia: Infinity,
    friction: Infinity,
    static: true
}));
world.add(new ji.RigidBody({
    pos: ji.V.V(width / 2, height),
    shape: new ji.shape.Box(0, 0, 0, width, 200),
    mass: Infinity,
    inertia: Infinity,
    friction: Infinity,
    static: true
}));
world.add(new ji.RigidBody({
    pos: ji.V.V(0, height / 2),
    shape: new ji.shape.Box(0, 0, 0, 200, height),
    mass: Infinity,
    inertia: Infinity,
    friction: Infinity,
    static: true
}));
world.add(new ji.RigidBody({
    pos: ji.V.V(width, height / 2),
    shape: new ji.shape.Box(0, 0, 0, 200, height),
    mass: Infinity,
    inertia: Infinity,
    friction: Infinity,
    static: true
}));

world.add(new ji.RigidBody({
    pos: ji.V.V(width / 2, height * 2 / 3 - 280),
    rot: 0,
    mass: 50,
    inertia: 50000,
    friction: 0.4,
    shape: new ji.shape.Box(0, 0, 0, 100, 100)
}));
world.add(new ji.RigidBody({
    pos: ji.V.V(width / 2, height * 2 / 3),
    mass: 100,
    inertia: 50000,
    friction: 0.8,
    shape: new ji.shape.Box(0, 0, 0, 100, 100)
}));

for(let i = 0; i < ropeLength; i += 1) {
    world.add(new ji.RigidBody({
        pos: ji.V.V(width / 3, 110 + i * 25),
        mass: 10,
        inertia: 100,
        friction: 0.4,
        shape: new ji.shape.Box(0, 0, 0, 10, 20)
    }));
}

world.add(new ji.RigidBody({
    pos: ji.V.V(130, height - 100 - 50 - 50 / 2 - 150 / 2),
    mass: 100,
    inertia: 100000,
    friction: 0.4,
    shape: new ji.shape.Box(0, 0, 0, 50, 150)
}));
world.add(new ji.RigidBody({
    pos: ji.V.V(130 + 200 / 2 - 50 / 2, height - 100 - 50),
    mass: 100,
    inertia: 100000,
    friction: 0.4,
    shape: new ji.shape.Box(0, 0, 0, 200, 50)
}));

for(let j = 0; j < 3; j += 1) {
    for(let k = 0; k < 3; k += 1) {
        const shape = [];
        for(let i = 0; i < 7; i += 1) {
            const a = i * Math.PI * 2 / 7 + Math.random() * 0.2 - 0.1;
            const r = Math.random() * 10 + 50;
            shape.push(ji.V.V(Math.cos(a) * r, Math.sin(a) * r));
        }
        world.add(new ji.RigidBody({
            pos: ji.V.V(width * 1/2 + 110 + j * 110, height / 2 + k * 110),
            rot: 0,
            mass: 100,
            inertia: 100000,
            friction: 1,
            shape: (j + k) % 2 === 0 ? new ji.shape.Convex(shape) : new ji.shape.Circle(0, 0, Math.random() * 10 + 50)
        }));
    }
}

let anchor = false, anchorBody;

let pt = Date.now();
function draw() {
    const t = Date.now();
    const dt = constrain(t - pt, 0, 20);
    background(255);
    const iterations = 6;
    for(let i = 0; i < iterations; i += 1) {
        world.simulate(dt * 0.001 / iterations);
    }
    world.clearConstraints();
    noFill(255);
    stroke(0);
    strokeWeight(1);
    if(mouseIsPressed && anchor) {
        if(mouseButton === LEFT) {
            world.addConstraint("anchor", {
                body: anchorBody,
                anchor,
                point: ji.V.V(mouseX, mouseY),
            });
        }
        else {
            world.addConstraint("anchor", {
                body: anchorBody,
                anchor,
                point: ji.V.V(mouseX, mouseY),
            }, "oscillation", 1, 1);
            const a = ji.V.mdot(anchorBody.transform(), anchor);
            line(mouseX, mouseY, a.x, a.y);
        }
        cursor("grabbing");
    }
    world.addConstraint("distance", {
        body1: world.bodies[4],
        anchor1: ji.V.V(0, 40),
        body2: world.bodies[5],
        anchor2: ji.V.V(0, -40),
        distance: 200
    })
    const a1 = ji.V.mdot(world.bodies[4].transform(), ji.V.V(0, 40));
    const a2 = ji.V.mdot(world.bodies[5].transform(), ji.V.V(0, -40));
    line(a1.x, a1.y, a2.x, a2.y);
    
    world.addConstraint("distance", {
        body1: world.bodies[0],
        anchor1: ji.V.V(-width / 2 + width / 3, 100),
        body2: world.bodies[6],
        anchor2: ji.V.V(0, -10),
        distance: 5
    });
    for(let i = 6; i < 6 + ropeLength - 1; i += 1) {
        world.addConstraint("distance", {
            body1: world.bodies[i],
            anchor1: ji.V.V(0, 10),
            body2: world.bodies[i + 1],
            anchor2: ji.V.V(0, -10),
            distance: 5
        });
        // const a1 = ji.V.mdot(world.bodies[i].transform(), ji.V.V(0, 10));
        // const a2 = ji.V.mdot(world.bodies[i + 1].transform(), ji.V.V(0, -10));
        // line(a1.x, a1.y, a2.x, a2.y);
    }
    world.addConstraint("union", {
        body1: world.bodies[6 + ropeLength],
        body2: world.bodies[6 + ropeLength + 1],
        pos: ji.V.V(200 / 2 - 50 / 2, 150 / 2 + 50 / 2),
        rot: 0
    });
    if(!mouseIsPressed) {
        cursor("default");
        anchor = false;
    }
    for(const body of world.bodies) {
        const shape = body.shape;
        noFill();
        if(!anchor && body.pointInside(ji.V.V(mouseX, mouseY))) {
            cursor("grab");
            if(mouseIsPressed) {
                anchor = ji.V.mdot(body.itransform(), ji.V.V(mouseX, mouseY));
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
            line(0, 0, shape.pos.x + shape.radius, shape.pos.y);
        }
        pop();
    }
    pt = t;
}