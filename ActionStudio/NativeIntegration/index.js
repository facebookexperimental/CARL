export async function initializeNativeIntegrationAsync() {
    return new Promise((resolve, reject) => {
        const script = document.createElement("script");
        script.type = "text/javascript";
        script.src = "/assets/carl.js";
        script.onload = () => CARL().then(carl => resolve(carl));
        document.body.appendChild(script);
    });
};
