var ji = (function () {

let DEBUG = false;
function debugMode(toggle = true) { DEBUG = toggle; }

const V = {
    V(x, y) {
        return { x, y };
    },
    sub(a, b) {
        return V.V(a.x - b.x, a.y - b.y);
    },
    isub(a, b) {
        a.x -= b.x;
        a.y -= b.y;
        return a;
    },
    add(a, b) {
        return V.V(a.x + b.x, a.y + b.y);
    },
    iadd(a, b) {
        a.x += b.x;
        a.y += b.y;
        return a;
    },
    inorm(v) {
        const d = 1 / Math.sqrt(v.x * v.x + v.y * v.y);
        v.x *= d;
        v.y *= d;
        return v;
    },
    scl(v, s) {
        return V.V(v.x * s, v.y * s);
    },
    iscl(v, s) {
        v.x *= s;
        v.y *= s;
        return v;
    },
    cross(a, b) {
        return (a.x * b.y) - (a.y * b.x);
    },
    dot(a, b) {
        return a.x * b.x + a.y * b.y;
    },
    mdot(m, v) {
        return V.V(v.x * m[0][0] + v.y * m[0][1] + m[0][2], v.x * m[1][0] + v.y * m[1][1] + m[1][2]);
        // return V.V(v.x * m[0].x + v.y * m[1].x + v.z * m[2].x, v.x * m[0].y + v.y * m[1].y + v.z * m[2].y, 1);
    },
    mag(v) {
        return Math.sqrt(v.x * v.x + v.y * v.y);
    },
    rot(v, theta) {
        const c = Math.cos(theta), s = Math.sin(theta);
        return V.V(v.x * c - v.y * s, v.x * s + v.y * c);
    },
    mrot(theta, t) {
        const c = Math.cos(theta), s = Math.sin(theta);
        return [
            [c, -s, t.x],
            [s,  c, t.y]
        ];
        // return [V.V(c, s, 0), V.V(-s, c, 0), V.V(t.x, t.y, 1)];
    },
};

class World {
    bodies = [];
    layers = {
        main: [],
    };
    constraints = [];
    constructor(gx = 0, gy = 0) {
        this.gravity = V.V(gx, gy);
    }
    add(body, layer = "main") {
        this.bodies.push(body);
        if(layer instanceof Array) {
            for(const l of layer) {
                this.layers[l] = this.bodies[l] || [];
                this.layers[l].push(body);
            }
        }
        else {
            this.layers[layer] = this.layers[layer] || [];
            this.layers[layer].push(body);
        }
    }
    /**
     * Adds a velocity constraint
     * 
     * `J * v + (beta/dt) * C + gamma * lambda = 0`  
     * `v += dt * M^-1 * J^T * lambda`
     * @param {"anchor" | "distance" | "union" | "joint"} type - "anchor," "distance," "union," or "joint"
     * @param {*} data - If "anchor," `{ body, anchor, point }`; if "distance," `{ body1, anchor1, body2, anchor2, distance }`; if "union," `{ body1, body2, pos, rot }`; if "joint," `{ body1, anchor1, body2, anchor2 }`
     * @param {"hard" | "oscillation" | "decay"} parameterization - "hard," "oscillation," or "decay"
     * @param {number} beta - If "hard," `beta`; if "oscillation," frequency; if "decay," half-life
     * @param {number} gamma - If "hard," `gamma * effective mass`; if "oscillation," half-life; if "decay," none
     */
    addConstraint(type, data, parameterization = ["hard", "oscillation", "decay"][0], beta = 0.3, gamma = 0) {
        this.constraints.push([type, data, parameterization, beta, gamma]);
    }
    collidable(body1, body2) {
        if(this.layers.main.includes(body1) || this.layers.main.includes(body2)) return true;
        for(const layer in this.layers) {
            if(this.layers[layer].includes(body1) && this.layers[layer].includes(body2)) return true;
        }
        return false;
    }
    remove(body) {
        for(const layer in this.layers) {
            if(this.layers[layer].includes(body)) this.layers[layer].splice(this.layers[layer].indexOf(body), 1);
        }
        this.bodies.splice(this.bodies.indexOf(body), 1);
    }
    solveConstraints(dt) {
        for(let i = 0; i < this.constraints.length; i += 1) {
            const constraint = this.constraints[i];
            let [type, data, parameterization, beta, gamma] = constraint;
            switch(parameterization) {
                case "hard": {
                    beta /= dt;
                    gamma = -1 / (1 + gamma / dt);
                } break;
                case "oscillation": {
                    const omega = 2 * Math.PI * beta, zeta = Math.log(2) / gamma;
                    const c = 2 * zeta, k = omega * omega;
                    beta = k / (c + dt * k);
                    gamma = -1 / (1 + 1 / (c * dt + k * dt * dt));
                } break;
                case "decay": {
                    const omega = Math.log(2) / beta, zeta = 1;
                    const c = 2 * omega * zeta, k = omega * omega;
                    beta = k / (c + dt * k);
                    gamma = -1 / (1 + 1 / (c * dt + k * dt * dt));
                } break;
            }
            switch(type) {
                case "anchor": {
                    let { body, anchor, point } = data;
                    anchor = V.rot(anchor, body.rot);
                    const og = [body.vel.x, body.vel.y, body.rvel];
                    const b = [body.pos.x + anchor.x - point.x, body.pos.y + anchor.y - point.y];
                    const J = [
                        [1, 0, -anchor.y],
                        [0, 1,  anchor.x]
                    ];
                    const iM = [
                        [body.imass,           0,            0],
                        [         0, body.imass,             0],
                        [         0,          0, body.iinertia]
                    ];
                    const iMJT = math.multiply(iM, math.transpose(J));
                    const solved = math.multiply(
                        iMJT,
                        math.multiply(
                            math.pinv(math.multiply(J, iMJT)),
                            math.multiply(
                                math.add(math.multiply(J, og), math.multiply(b, beta)),
                                gamma
                            )
                        )
                    );
                    body.vel.x += solved[0];
                    body.vel.y += solved[1];
                    body.rvel += solved[2];
                } break;
                case "distance": {
                    let { body1, anchor1, body2, anchor2, distance = 0 } = data;
                    anchor1 = V.mdot(body1.transform(), anchor1);
                    anchor2 = V.mdot(body2.transform(), anchor2);
                    let dx = anchor2.x - anchor1.x, dy = anchor2.y - anchor1.y;
                    const og = [body1.vel.x, body1.vel.y, body1.rvel, body2.vel.x, body2.vel.y, body2.rvel];
                    let d = Math.sqrt(dx * dx + dy * dy);
                    const b = d - distance;
                    d = 1 / d;
                    dx *= d;
                    dy *= d;
                    const J = [[
                        -dx,
                        -dy,
                        dx * (anchor1.y - body1.pos.y) - dy * (anchor1.x - body1.pos.x),

                        dx,
                        dy,
                        dy * (anchor2.x - body2.pos.x) - dx * (anchor2.y - body2.pos.y)
                    ]];
                    
                    const iM = [
                        [body1.imass, 0, 0, 0, 0, 0],
                        [0, body1.imass, 0, 0, 0, 0],
                        [0, 0, body1.iinertia, 0, 0, 0],
                        [0, 0, 0, body2.imass, 0, 0],
                        [0, 0, 0, 0, body2.imass, 0],
                        [0, 0, 0, 0, 0, body2.iinertia]
                    ];
                    const iMJT = math.multiply(iM, math.transpose(J));
                    const solved = math.multiply(
                        iMJT,
                        math.multiply(
                            math.pinv(math.multiply(J, iMJT)),
                            math.multiply(
                                math.add(math.multiply(J, og), math.multiply(b, beta)),
                                gamma
                            )
                        )
                    );
                    body1.vel.x += solved[0];
                    body1.vel.y += solved[1];
                    body1.rvel += solved[2];
                    body2.vel.x += solved[3];
                    body2.vel.y += solved[4];
                    body2.rvel += solved[5];
                } break;
                case "union": {
                    let { body1, body2, pos, rot } = data;
                    const cb1 = Math.cos(body1.rot), sb1 = Math.sin(body1.rot);
                    const og = [body1.vel.x, body1.vel.y, body1.rvel, body2.vel.x, body2.vel.y, body2.rvel];
                    const b = [
                        body1.pos.x - body2.pos.x + pos.x * cb1 - pos.y * sb1,
                        body1.pos.y - body2.pos.y + pos.x * sb1 + pos.y * cb1,
                        body2.rot - body1.rot - rot
                    ];
                    const J = [
                        [1, 0, -pos.x * sb1 - pos.y * cb1, -1, 0, 0],
                        [0, 1, pos.x * cb1 - pos.y * sb1, 0, -1, 0],
                        [0, 0, -1, 0, 0, 1]
                    ];
                    const iM = [
                        [body1.imass, 0, 0, 0, 0, 0],
                        [0, body1.imass, 0, 0, 0, 0],
                        [0, 0, body1.iinertia, 0, 0, 0],
                        [0, 0, 0, body2.imass, 0, 0],
                        [0, 0, 0, 0, body2.imass, 0],
                        [0, 0, 0, 0, 0, body2.iinertia]
                    ];
                    const iMJT = math.multiply(iM, math.transpose(J));
                    const solved = math.multiply(
                        iMJT,
                        math.multiply(
                            math.pinv(math.multiply(J, iMJT)),
                            math.multiply(
                                math.add(math.multiply(J, og), math.multiply(b, beta)),
                                gamma
                            )
                        )
                    );
                    body1.vel.x += solved[0];
                    body1.vel.y += solved[1];
                    body1.rvel += solved[2];
                    body2.vel.x += solved[3];
                    body2.vel.y += solved[4];
                    body2.rvel += solved[5];
                } break;
                case "joint": {
                    let { body1, anchor1, body2, anchor2 } = data;
                    anchor1 = V.mdot(body1.transform(), anchor1);
                    anchor2 = V.mdot(body2.transform(), anchor2);
                    const og = [body1.vel.x, body1.vel.y, body1.rvel, body2.vel.x, body2.vel.y, body2.rvel];
                    const b = [
                        anchor1.x - anchor2.x,
                        anchor1.y - anchor2.y
                    ];
                    const J = [
                        [1, 0, -(anchor1.y - body1.pos.y), -1, 0, (anchor2.y - body2.pos.y)],
                        [0, 1, (anchor1.x - body1.pos.x), 0, -1, -(anchor2.x - body2.pos.x)]
                    ];
                    const iM = [
                        [body1.imass, 0, 0, 0, 0, 0],
                        [0, body1.imass, 0, 0, 0, 0],
                        [0, 0, body1.iinertia, 0, 0, 0],
                        [0, 0, 0, body2.imass, 0, 0],
                        [0, 0, 0, 0, body2.imass, 0],
                        [0, 0, 0, 0, 0, body2.iinertia]
                    ];
                    const iMJT = math.multiply(iM, math.transpose(J));
                    const solved = math.multiply(
                        iMJT,
                        math.multiply(
                            math.pinv(math.multiply(J, iMJT)),
                            math.multiply(
                                math.add(math.multiply(J, og), math.multiply(b, beta)),
                                gamma
                            )
                        )
                    );
                    body1.vel.x += solved[0];
                    body1.vel.y += solved[1];
                    body1.rvel += solved[2];
                    body2.vel.x += solved[3];
                    body2.vel.y += solved[4];
                    body2.rvel += solved[5];
                } break;
            }
        }
    }
    simulate(dt, contactBias = 0.1) {
        // integrate velocity
        for(let i = 0; i < this.bodies.length; i += 1) {
            const body = this.bodies[i];

            if(!body.static) V.iadd(body.vel, V.scl(this.gravity, dt));
        }

        // velocity constraints
        for(let i = 0; i < this.bodies.length - 1; i += 1) {
            const body = this.bodies[i];
            for(let j = i + 1; j < this.bodies.length; j += 1) {
                let body2 = this.bodies[j];
                let body1 = body.static ? body2 : body;
                if(!this.collidable(body1, body2)) continue;
                body2 = body.static ? body : body2;
                if(body2.static && body1.static) continue;
                const manifolds = getManifolds(body1, body2);
                if(manifolds) {
                    const { contactPoints, referenceNormal, incidentNormal, referenceBody, incidentBody } = manifolds[0];

                    // normal
                    let og = [referenceBody.vel.x, referenceBody.vel.y, referenceBody.rvel, incidentBody.vel.x, incidentBody.vel.y, incidentBody.rvel];
                    let Jn = [];
                    let bn = [];
                    let collidingPoints = [];
                    for(const contact of contactPoints) {
                        const { pointOnIncident, pointOnReference, penetration } = contact;
                        if(DEBUG) {
                            line(pointOnReference.x, pointOnReference.y, referenceBody.pos.x, referenceBody.pos.y);
                            line(pointOnIncident.x, pointOnIncident.y, incidentBody.pos.x, incidentBody.pos.y);
                        }

                        const restitution = Math.max(body1.restitution, body2.restitution);
                        let Jni = [
                            -referenceNormal.x,
                            -referenceNormal.y,
                            referenceNormal.x * (pointOnReference.y - referenceBody.pos.y) - referenceNormal.y * (pointOnReference.x - referenceBody.pos.x),
                            
                            referenceNormal.x,
                            referenceNormal.y,
                            referenceNormal.y * (pointOnIncident.x - incidentBody.pos.x) - referenceNormal.x * (pointOnIncident.y - incidentBody.pos.y)
                        ],
                            bni = (penetration * contactBias / dt + restitution * math.dot(Jni, og)) * (referenceBody === referenceBody ? -1 : 1);
                        if(math.dot(Jni, og) < -bni) {
                            Jn.push(Jni);
                            bn.push(bni);
                            collidingPoints.push(contact);
                        }
                    }
                    let lambdaN = 0;
                    if(Jn.length) {
                        const iM = [
                            [referenceBody.imass, 0, 0, 0, 0, 0],
                            [0, referenceBody.imass, 0, 0, 0, 0],
                            [0, 0, referenceBody.iinertia, 0, 0, 0],
                            [0, 0, 0, incidentBody.imass, 0, 0],
                            [0, 0, 0, 0, incidentBody.imass, 0],
                            [0, 0, 0, 0, 0, incidentBody.iinertia]
                        ];
                        const iMJT = math.multiply(iM, math.transpose(Jn));
                        lambdaN = math.multiply(math.pinv(math.multiply(Jn, iMJT)), math.multiply(math.add(math.multiply(Jn, og), bn), -1));
                        const solved = math.multiply(iMJT, lambdaN);
                        referenceBody.vel.x += solved[0];
                        referenceBody.vel.y += solved[1];
                        referenceBody.rvel += solved[2];
                        incidentBody.vel.x += solved[3];
                        incidentBody.vel.y += solved[4];
                        incidentBody.rvel += solved[5];
                    }

                    // friction
                    if(lambdaN) {
                        let og = [referenceBody.vel.x, referenceBody.vel.y, referenceBody.rvel, incidentBody.vel.x, incidentBody.vel.y, incidentBody.rvel];
                        let Jf = [];
                        const tangent = V.V(-referenceNormal.y, referenceNormal.x);
                        const friction = Math.min(body1.friction, body2.friction);
                        for(const contact of collidingPoints) {
                            const { pointOnIncident, pointOnReference, penetration } = contact;

                            let Jfi = [
                                -tangent.x,
                                -tangent.y,
                                tangent.x * (pointOnReference.y - referenceBody.pos.y) - tangent.y * (pointOnReference.x - referenceBody.pos.x),
                                
                                tangent.x,
                                tangent.y,
                                tangent.y * (pointOnIncident.x - incidentBody.pos.x) - tangent.x * (pointOnIncident.y - incidentBody.pos.y)
                            ];
                            Jf.push(Jfi);
                        }
                        if(Jf.length) {
                            const iM = [
                                [referenceBody.imass, 0, 0, 0, 0, 0],
                                [0, referenceBody.imass, 0, 0, 0, 0],
                                [0, 0, referenceBody.iinertia, 0, 0, 0],
                                [0, 0, 0, incidentBody.imass, 0, 0],
                                [0, 0, 0, 0, incidentBody.imass, 0],
                                [0, 0, 0, 0, 0, incidentBody.iinertia]
                            ];
                            const iMJT = math.multiply(iM, math.transpose(Jf));
                            let lambdaF = math.multiply(math.pinv(math.multiply(Jf, iMJT)), math.multiply(math.multiply(Jf, og), -1));
                            lambdaF = lambdaF.map((a, i) => Math.max(Math.min(a, friction * Math.abs(lambdaN[0])), -friction * Math.abs(lambdaN[0])));
                            const solved = math.multiply(iMJT, lambdaF);
                            referenceBody.vel.x += solved[0];
                            referenceBody.vel.y += solved[1];
                            referenceBody.rvel += solved[2];
                            incidentBody.vel.x += solved[3];
                            incidentBody.vel.y += solved[4];
                            incidentBody.rvel += solved[5];
                        }
                    }
                }
            }
        }
        this.solveConstraints(dt);

        // integrate position
        for(let i = 0; i < this.bodies.length; i += 1) {
            const body = this.bodies[i];
            if(!body.static) {
                V.iadd(body.pos, V.scl(body.vel, dt));
                body.rot += body.rvel * dt;
            }
        }
        // position constraints
    }
    clearConstraints() {
        this.constraints.length = 0;
    }
}

class RigidBody {
    #mass; #imass; #inertia; #iinertia;
    vel = V.V(0, 0);
    rvel = 0;
    /**
     * Creates a new RigidBody
     * @param {number} x - x coordinate of center of mass
     * @param {number} y - y coordinate of CoM
     * @param {*} shape - Physical shape
     * @param {number} mass - Mass
     * @param {number} rot - Rotation
     * @param {number} inertia - Moment of inertia
     * @param {boolean} s - Static
     */
    constructor({ pos, shape, mass = 1, rot = 0, inertia = 1, restitution = 0, friction = 0, static: s = false }) {
        this.pos = pos;
        this.rot = rot;
        this.mass = mass;
        this.inertia = inertia;
        this.shape = shape;
        this.restitution = restitution;
        this.friction = friction;
        this.static = s;
    }

    /**
     * @param {number} mass
     */
    set mass(mass) {
        this.#mass = mass;
        this.#imass = 1 / mass;
    }
    get mass() { return this.#mass }
    /**
     * @param {number} imass
     */
    set imass(imass) {
        this.#imass = imass;
        this.#mass = 1 / imass;
    }
    get imass() { return this.#imass; }
    /**
     * @param {number} inertia
     */
    set inertia(inertia) {
        this.#inertia = inertia;
        this.#iinertia = 1 / inertia;
    }
    get inertia() { return this.#inertia }
    /**
     * @param {number} iinertia
     */
    set iinertia(iinertia) {
        this.#iinertia = iinertia;
        this.#inertia = 1 / iinertia;
    }
    get iinertia() { return this.#iinertia; }

    transform() { // local to world
        return V.mrot(this.rot, this.pos);
    }
    itransform() { // world to local
        return V.mrot(-this.rot, V.rot(this.pos, Math.PI - this.rot));
    }
    pointInside(worldPoint) {
        return this.shape.pointInside(V.mdot(this.itransform(), worldPoint));
    }
    velAtPoint(worldPoint) {
        return V.iadd(V.iscl(V.V(this.pos.y - worldPoint.y, worldPoint.x - this.pos.x), this.rvel), this.vel);
    }
}

class Convex {
    /**
     * Creates a convex polygon out of an array of vertices. The original array will be modified to be clockwise (up positive).
     * @param {*} verts 
     */
    constructor(verts) {
        this.verts = verts;
        for(let i = 0; i < verts.length; i += 1) {
            const c = V.cross(V.sub(verts[i + 1], verts[i]), V.sub(verts[i + 2], verts[i + 1]));
            if(c < 0) return;
            if(c > 0) {
                verts.reverse();
                return;
            }
        }
        this.bounding = verts[0].x * verts[0].x + verts[0].y * verts[0].y;
        for(let i = 1; i < verts.length; i += 1) {
            const d = verts[i].x * verts[i].x + verts[i].y * verts[i].y;
            if(d < this.bounding) this.bounding = d;
        }
        this.bounding = Math.sqrt(this.bounding);
    }
    world(transform) {
        const verts = [];
        for(let i = 0; i < this.verts.length; i += 1) {
            verts.push(V.mdot(transform, this.verts[i]));
        }
        return verts;
    }
    pointInside(localPoint) {
        for(let i = 0; i < this.verts.length; i += 1) {
            const vert1 = this.verts[i], vert2 = this.verts[(i + 1) % this.verts.length];
            const normal = V.V(vert1.y - vert2.y, vert2.x - vert1.x);
            if(V.dot(V.sub(localPoint, vert1), normal) > 0) return false;
        }
        return true;
    }
}
class Box extends Convex {
    constructor(x = 0, y = 0, rot = 0, w = 1, h = 1) {
        const pos = V.V(x, y);
        const m = V.mrot(rot, pos);
        super([
            V.mdot(m, V.V(-w / 2, h / 2)),
            V.mdot(m, V.V(w / 2, h / 2)),
            V.mdot(m, V.V(w / 2, -h / 2)),
            V.mdot(m, V.V(-w / 2, -h / 2))
        ]);
        this.pos = pos;
        this.size = V.V(w, h);
    }
}

class Circle {
    constructor(x = 0, y = 0, radius = 10) {
        this.pos = V.V(x, y);
        this.radius = radius;
        this.bounding = radius + V.mag(this.pos);
    }
    pointInside(localPoint) {
        return localPoint.x * localPoint.x + localPoint.y * localPoint.y < this.radius * this.radius;
    }
    world(transform) {
        return {
            pos: V.mdot(transform, this.pos),
            radius: this.radius
        };
    }
}

/**
 * Gets contact manifolds between two bodies
 * @param {RigidBody} body1 
 * @param {RigidBody} body2 
 * @returns {[[[{ pointOnIncident, pointOnReference, penetration }], referenceNormal, incidentNormal, referenceBody, incidentBody]]} - 
 */
function getManifolds(body1, body2) {
    if((body2.pos.x - body1.pos.x) * (body2.pos.x - body1.pos.x) + (body2.pos.y - body1.pos.y) * (body2.pos.y - body1.pos.y) > (body1.shape.bounding + body2.shape.bounding) * (body1.shape.bounding + body2.shape.bounding)) return false;
    let shapesSwapped = 0;
    let shape1 = body1.shape.world(body1.transform()),
        shape2 = body2.shape.world(body2.transform());
    if(body2.shape instanceof Convex) {
        shapesSwapped = 1;
        [shape1, shape2] = [shape2, shape1];
        [body1, body2] = [body2, body1];
    }
    if(body1.shape instanceof Circle) {
        if(body2.shape instanceof Circle) {
            if((shape2.pos.x - shape1.pos.x) * (shape2.pos.x - shape1.pos.x) + (shape2.pos.y - shape1.pos.y) * (shape2.pos.y - shape1.pos.y) > (shape1.radius + shape2.radius) * (shape1.radius * shape2.radius)) return false;
            const v = V.sub(body1.pos, body2.pos);
            let d = V.mag(v);
            const penetration = shape1.radius + shape2.radius - d, referenceBody = body1, incidentBody = body2;
            const referenceNormal = V.iscl(V.sub(incidentBody.pos, referenceBody.pos), 1 / d);
            const contactPoints = [{
                pointOnIncident: V.add(incidentBody.pos, V.scl(referenceNormal, -incidentBody.shape.radius)),
                pointOnReference: V.add(referenceBody.pos, V.scl(referenceNormal, referenceBody.shape.radius)),
                penetration: penetration
            }];
            return [{ contactPoints, referenceNormal, incidentNormal: V.scl(referenceNormal, -1), referenceBody, incidentBody }];
        }
    }
    if(body1.shape instanceof Convex) {
        if(body2.shape instanceof Circle) {
            let minPenetration = Infinity, referenceBody = body1, incidentBody = body2, referenceNormal = 0, referenceSupport = 0;
            let closest = Infinity;
            let centerOutside = -1;
            for(let i = 0; i < shape1.length; i += 1) {
                const vert1 = shape1[i], vert2 = shape1[(i + 1) % shape1.length];
                const vector = V.sub(shape2.pos, vert1);
                const edge = V.sub(vert2, vert1);
                const d = V.mag(edge), id = 1 / d;
                const normal = V.V(-edge.y * id, edge.x * id);
                const t = V.dot(vector, edge) * id;
                const penetration = shape2.radius - V.dot(normal, vector);
                if(penetration < 0) return false;
                if(penetration < shape2.radius) centerOutside = 1;
                const projected = V.add(vert1, V.scl(edge, id * Math.max(0, Math.min(d, t))));
                const distance = (projected.x - shape2.pos.x) * (projected.x - shape2.pos.x) + (projected.y - shape2.pos.y) * (projected.y - shape2.pos.y);
                if(distance < closest) {
                    closest = distance;
                    referenceSupport = projected;
                }
            }
            referenceNormal = V.inorm(V.iscl(V.sub(shape2.pos, referenceSupport), centerOutside));
            minPenetration = shape2.radius - V.dot(referenceNormal, V.sub(shape2.pos, referenceSupport));
            const contactPoints = [{
                pointOnIncident: V.add(referenceSupport, V.scl(referenceNormal, -minPenetration)),
                pointOnReference: referenceSupport,
                penetration: minPenetration
            }];
            return [{ contactPoints, referenceNormal, incidentNormal: V.scl(referenceNormal, -1), referenceBody, incidentBody }];
        }
        else if(body2.shape instanceof Convex) {
            // get reference face and penetration
            let minPenetration = Infinity, referenceBody = body1, referenceFace = 0, referenceNormal = 0, incidentSupport = 0;
            for(let i = 0; i < shape1.length; i += 1) {
                const vert1 = shape1[i], vert2 = shape1[(i + 1) % shape1.length];
                const normal = V.inorm(V.V(vert1.y - vert2.y, vert2.x - vert1.x));
                let support = shape2[0], penetration = V.dot(V.sub(vert1, support), normal);
                for(let j = 1; j < shape2.length; j += 1) {
                    const newPenetration = V.dot(V.sub(vert1, shape2[j]), normal);
                    if(newPenetration > penetration) {
                        penetration = newPenetration;
                        support = shape2[j];
                    }
                }
                if(penetration < 0) return false;
                if(penetration < minPenetration) {
                    minPenetration = penetration;
                    referenceBody = body1;
                    incidentBody = body2;
                    referenceFace = i;
                    incidentSupport = support;
                    referenceNormal = normal;
                }
            }
            for(let i = 0; i < shape2.length; i += 1) {
                const vert1 = shape2[i], vert2 = shape2[(i + 1) % shape2.length];
                const normal = V.inorm(V.V(vert1.y - vert2.y, vert2.x - vert1.x));
                let support = shape1[0], penetration = V.dot(V.sub(vert1, support), normal);
                for(let j = 1; j < shape1.length; j += 1) {
                    const newPenetration = V.dot(V.sub(vert1, shape1[j]), normal);
                    if(newPenetration > penetration) {
                        penetration = newPenetration;
                        support = shape1[j];
                    }
                }
                if(penetration < 0) return false;
                if(penetration < minPenetration) {
                    minPenetration = penetration;
                    referenceBody = body2;
                    incidentBody = body1;
                    referenceFace = i;
                    incidentSupport = support;
                    referenceNormal = normal;
                }
            }
            
            // get incident face
            if(referenceBody === body2) {
                shapesSwapped = 1 - shapesSwapped;
                [shape1, shape2] = [shape2, shape1];
                [body1, body2] = [body2, body1];
            }
            let incidentFace = 0, incidentNormal = 0, minDot = Infinity;
            for(let i = 0; i < shape2.length; i += 1) {
                const vert1 = shape2[i], vert2 = shape2[(i + 1) % shape2.length];
                const normal = V.inorm(V.V(vert1.y - vert2.y, vert2.x - vert1.x));
                const dot = V.dot(normal, referenceNormal);
                if(dot < minDot) {
                    minDot = dot;
                    incidentFace = i;
                    incidentNormal = normal;
                }
            }

            // clip incident face to sides of reference face
            let inc1 = shape2[incidentFace], inc2 = shape2[(incidentFace + 1) % shape2.length];
            let ref1 = shape1[referenceFace], ref2 = shape1[(referenceFace + 1) % shape1.length];
            let denom = 1 / (referenceNormal.x * (inc1.y - inc2.y) + referenceNormal.y * (inc2.x - inc1.x));
            let x1 = Math.max(
                (referenceNormal.x * (ref1.y - inc2.y) + referenceNormal.y * (inc2.x - ref1.x)) * denom,
                (referenceNormal.x * (ref2.y - inc2.y) + referenceNormal.y * (inc2.x - ref2.x)) * denom
            );
            inc1 = V.add(inc2, V.iscl(V.sub(inc1, inc2), Math.min(1, x1 > 0 ? x1 : 1)));

            denom = 1 / (referenceNormal.x * (inc2.y - inc1.y) + referenceNormal.y * (inc1.x - inc2.x));
            x1 = Math.max(
                (referenceNormal.x * (ref1.y - inc1.y) + referenceNormal.y * (inc1.x - ref1.x)) * denom,
                (referenceNormal.x * (ref2.y - inc1.y) + referenceNormal.y * (inc1.x - ref2.x)) * denom
            );
            inc2 = V.add(inc1, V.iscl(V.sub(inc2, inc1), Math.min(1, x1 > 0 ? x1 : 1)));

            const contactPoints = [];
            let dot = V.dot(referenceNormal, V.sub(inc1, ref1));
            if(dot <= 0) contactPoints.push({
                pointOnIncident: V.V(inc1.x, inc1.y),
                pointOnReference: V.V(inc1.x - referenceNormal.x * dot, inc1.y - referenceNormal.y * dot),
                penetration: -dot
            });
            dot = V.dot(referenceNormal, V.sub(inc2, ref1));
            if(dot <= 0) contactPoints.push({
                pointOnIncident: V.V(inc2.x, inc2.y),
                pointOnReference: V.V(inc2.x - referenceNormal.x * dot, inc2.y - referenceNormal.y * dot),
                penetration: -dot
            });

            if(DEBUG) {
                strokeWeight(2);
                line(incidentSupport.x, incidentSupport.y, incidentSupport.x + referenceNormal.x * minPenetration, incidentSupport.y + referenceNormal.y * minPenetration);
                let vert = V.iscl(V.add(shape2[incidentFace], shape2[(incidentFace + 1) % shape2.length]), 1/2);
                stroke(255, 0, 0);
                line(vert.x, vert.y, vert.x + incidentNormal.x * 50, vert.y + incidentNormal.y * 50);
                vert = V.iscl(V.add(shape1[referenceFace], shape1[(referenceFace + 1) % shape1.length]), 1/2);
                stroke(0);
                line(vert.x, vert.y, vert.x + referenceNormal.x * 50, vert.y + referenceNormal.y * 50);
                stroke(255, 0, 0);
                strokeWeight(5);
                line(inc1.x, inc1.y, inc2.x, inc2.y);
                strokeWeight(10);
                stroke(0, 255, 0);
                point(inc1.x, inc1.y);
                stroke(0, 0, 255);
                point(inc2.x, inc2.y);
                for(const p of contactPoints) {
                    strokeWeight(1);
                    line(p.pointOnIncident.x, p.pointOnIncident.y, p.pointOnReference.x, p.pointOnReference.y);
                }
            }
            return [{ contactPoints, referenceNormal, incidentNormal, referenceBody, incidentBody }];
        }
    }
}

return {
    World, RigidBody,
    shape: {
        Box, Convex, Circle
    },
    debugMode,
    getManifolds,
    V,
}

})();