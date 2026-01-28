export function testNativeIntegration() {
    const script = document.createElement("script");
    script.type = "text/javascript";
    script.src = "assets/carl.js";
    script.onload = () => {
        console.log("CARL WASM loaded! In theory!");
    };
    document.body.appendChild(script);

    console.log("Native integration test function called.");
};
