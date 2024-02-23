import * as THREE from 'three';
import Stats from 'three/addons/libs/stats.module.js';
import { GUI } from 'three/addons/libs/lil-gui.module.min.js';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { OBJLoader } from 'three/addons/loaders/OBJLoader.js';

// Assuming container IDs are like 'container1', 'container2', etc.
const objFiles = [
    '../../assets/output/quadratic_error_simplification/armadillo_20perc.obj',
    '../../assets/output/quadratic_error_simplification/armadillo_40perc.obj',
	'../../assets/output/quadratic_error_simplification/armadillo_60perc.obj',
	'../../assets/output/quadratic_error_simplification/armadillo_90perc.obj',
	'../../assets/output/quadratic_error_simplification/bunny_20perc.obj',
	'../../assets/output/quadratic_error_simplification/bunny_40perc.obj',
	'../../assets/output/quadratic_error_simplification/bunny_60perc.obj',
	'../../assets/output/quadratic_error_simplification/bunny_90perc.obj',
	'../../assets/output/quadratic_error_simplification/cow_20perc.obj',
	'../../assets/output/quadratic_error_simplification/cow_40perc.obj',
	'../../assets/output/quadratic_error_simplification/cow_60perc.obj',
	'../../assets/output/quadratic_error_simplification/cow_90perc.obj',
	'../../assets/output/loop_subdivision/cube_subdivided.obj',

    // Add more file paths as needed
];

let visualizations = [];

function initVisualization(containerId, objFilePath) {
    let container = document.getElementById(containerId);
    container.style.position = 'relative';

    let scene = new THREE.Scene();
    scene.background = new THREE.Color(0xffffff);
    let camera = new THREE.PerspectiveCamera(75, window.innerWidth / (window.innerHeight * 0.5), 0.1, 1000);
    
    let renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight * 0.5);
    container.appendChild(renderer.domElement);

    let controls = new OrbitControls(camera, renderer.domElement);
    controls.minDistance = 2;
    controls.maxDistance = 10;
    
    let dirlight = new THREE.DirectionalLight(0xffffff, 0.5);
    dirlight.position.set(0, 0, 1);
    scene.add(dirlight);

    let ambientLight = new THREE.AmbientLight(0x404040, 2);
    scene.add(ambientLight);

    let loader = new OBJLoader();
    loader.load(
        objFilePath,
        function (object) {
            let mesh = object.children[0];
            mesh.material = new THREE.MeshPhongMaterial({color: 0x999999});
            mesh.position.set(0, 0, 0);
            scene.add(mesh);

            // Initialize GUI after the object is loaded
            initGUI(container, mesh);
        },
        function (xhr) {
            console.log((xhr.loaded / xhr.total * 100) + '% loaded');
        },
        function (error) {
            console.log('An error happened: ' + error);
        }
    );

    camera.position.z = 5;

    visualizations.push({container, scene, camera, renderer, controls});
}

function initGUI(container, mesh) {
    let gui = new GUI({ autoPlace: false });
    gui.add(mesh.position, 'x', -1, 1);
    gui.add(mesh.position, 'y', -1, 1);
    gui.add(mesh.position, 'z', -1, 1);
    container.appendChild(gui.domElement);
    gui.domElement.style.position = 'absolute';
    gui.domElement.style.top = '0px';
    gui.domElement.style.right = '0px';
}

function animate() {
    requestAnimationFrame(animate);
    visualizations.forEach(vis => {
        vis.renderer.render(vis.scene, vis.camera);
    });
}

function onWindowResize() {
    visualizations.forEach(vis => {
        vis.camera.aspect = window.innerWidth / (window.innerHeight * 0.5);
        vis.camera.updateProjectionMatrix();
        vis.renderer.setSize(window.innerWidth, window.innerHeight * 0.5);
    });
}

window.addEventListener('resize', onWindowResize, false);

// Initialize each visualization
objFiles.forEach((file, index) => initVisualization(`container${index + 1}`, file));

animate();
